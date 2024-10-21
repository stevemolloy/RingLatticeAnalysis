#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "lib.h"
#include "elements.h"

#define M_PI (acos(-1.0))
#define ANGLE_EPSILON (0.1*M_PI/180)

#include "stb_ds.h"

// TODO: Rationalise memory management. In some places it is dynamic arrays, in other an arena is used. Time to straighten this out.
// TODO: Spin off a few more functions for the library.  The "main" function should be a list of calls to library functions.
// TODO: Rationalise organisation of the libraries.
// TODO: Handle off-energy dynamics.

void usage(const char* program) {
  fprintf(stderr, "USAGE: %s [-p <periodicity>] [-h <harmonic number>] [-E <kinetic energy>] lattice_file\n", program);
}

bool lattice_is_closed(double total_angle) {
  return !(fabs(total_angle - 2*M_PI) > ANGLE_EPSILON);
}

int main(int argc, char **argv) {
  CommandLineArgs args = get_clargs(argc, argv);

  Element *line = {0};
  generate_lattice(args.file_path, &line);

  double line_length = calculate_line_length(line);
  double total_length = line_length * args.periodicity;

  double line_angle = calculate_line_angle(line);
  double total_angle = line_angle * args.periodicity;

  double line_matrix[BEAM_DOFS*BEAM_DOFS] = {0};
  get_line_matrix(line_matrix, line);

  printf("\nSummary of the lattice defined in %s\n\n", args.file_path);

  bool closed_system = lattice_is_closed(total_angle);

  if (!closed_system) {
    printf("Total bending angle (%0.3f degrees) is not within %0.1e degrees of a circle.\n", 
           radians_to_degrees(total_angle), radians_to_degrees(ANGLE_EPSILON));
    printf("System does not close, so not calculating ring parameters.\n\n");
  }

  printf("Periodicity: %zu\n", args.periodicity);
  printf("Harmonic number: %d\n", args.harmonic_number);
  printf("Number of elements in the line: %td\n", arrlen(line));
  if (args.periodicity != 1) {
    printf("Total length of the lattice: %0.3f m (%0.3f m for the line)\n", total_length, line_length);
  } else {
    printf("Total length of the line: %f m\n", line_length);
  }
  if (args.periodicity != 1) {
    printf("Total bending angle of the lattice: %0.3f deg (%0.3f deg for the line))\n", 
           radians_to_degrees(total_angle), radians_to_degrees(line_angle));
  } else {
    printf("Total bending angle of the line: %0.3f degrees\n", radians_to_degrees(line_angle));
  }

  if (closed_system) {
    double rf_freq = C / (total_length / args.harmonic_number);
    printf("RF frequency for a hamonic number of %d: %0.6f MHz\n", args.harmonic_number, rf_freq/1e6);
  }

  printf("\nTotal matrix, R, for the line is:\n");
  rmatrix_print(stdout, line_matrix);

  double total_matrix[BEAM_DOFS*BEAM_DOFS] = {0};
  apply_matrix_n_times(total_matrix, line_matrix, args.periodicity);
  printf("\nTotal matrix, R, for the full system:\n");
  rmatrix_print(stdout, total_matrix);

  if (closed_system) {
    double x_trace = (total_matrix[0*BEAM_DOFS + 0] + total_matrix[1*BEAM_DOFS + 1]);
    double y_trace = (total_matrix[2*BEAM_DOFS + 2] + total_matrix[3*BEAM_DOFS + 3]);

    const double R11 = total_matrix[0*BEAM_DOFS + 0];
    const double R12 = total_matrix[0*BEAM_DOFS + 1];
    const double R22 = total_matrix[1*BEAM_DOFS + 1];
    const double R33 = total_matrix[2*BEAM_DOFS + 2];
    const double R34 = total_matrix[2*BEAM_DOFS + 3];
    const double R44 = total_matrix[3*BEAM_DOFS + 3];
    const double R16 = total_matrix[0*BEAM_DOFS + 5];
    const double R56 = total_matrix[4*BEAM_DOFS + 5];

    const double eta_x  = R16 / (1 - R11);
    const double phi_x = acos((R11 + R22) / 2.0f);
    const double phi_y = acos((R33 + R44) / 2.0f);
    const double beta_x = R12 / sin(phi_x);
    const double beta_y = R34 / sin(phi_y);

    double *Ss = {0};
    double *element_etas = {0};
    double *element_etaps = {0};
    double *element_beta_xs = {0};
    double *element_beta_ys = {0};
    double *element_curlyH = {0};

    double S = 0.0;
    double eta_vec[3] = {eta_x, 0.0f, 1.0f};
    double twiss_x_vec[3] = {beta_x, 0.0f, 1/beta_x};
    double twiss_y_vec[3] = {beta_y, 0.0f, 1/beta_y};

    arrput(Ss, 0.0);
    arrput(element_etas, eta_vec[0]);
    arrput(element_etaps, eta_vec[1]);
    arrput(element_beta_xs, beta_x);
    arrput(element_beta_ys, beta_y);
    arrput(element_curlyH, get_curlyH(eta_vec[0], eta_vec[1], twiss_x_vec[0], twiss_x_vec[1]));
    for (size_t i=0; i<arrlenu(line); i++) {
      S += element_length(line[i]);
      arrput(Ss, S);

      double temp_twiss_x_vec[3] = {0};
      matrix_multiply(line[i].twiss_prop_matrix_x, twiss_x_vec, temp_twiss_x_vec, 3, 3, 3, 1);
      memcpy(twiss_x_vec, temp_twiss_x_vec, 3*sizeof(double));
      arrput(element_beta_xs, twiss_x_vec[0]);

      double temp_twiss_y_vec[3] = {0};
      matrix_multiply(line[i].twiss_prop_matrix_y, twiss_y_vec, temp_twiss_y_vec, 3, 3, 3, 1);
      memcpy(twiss_y_vec, temp_twiss_y_vec, 3*sizeof(double));
      arrput(element_beta_ys, twiss_y_vec[0]);

      double temp_eta_vec[3] = {0};
      matrix_multiply(line[i].eta_prop_matrix, eta_vec, temp_eta_vec, 3, 3, 3, 1);
      memcpy(eta_vec, temp_eta_vec, 3*sizeof(double));
      arrput(element_etas, eta_vec[0]);
      arrput(element_etaps, eta_vec[1]);

      arrput(element_curlyH, get_curlyH(eta_vec[0], eta_vec[1], twiss_x_vec[0], twiss_x_vec[1]));
    }

    if (args.save_twiss) {
      FILE *twiss_file = fopen(args.twiss_filename, "w");
      fprintf(twiss_file, "S / m, beta_x / m, beta_y / m, eta_x / m\n");
      for (size_t i=0; i<arrlenu(line); i++) {
        fprintf(twiss_file, "%0.6e, %0.6e, %0.6e, %0.6e\n", Ss[i], element_beta_xs[i], element_beta_ys[i], element_etas[i]);
      }
      fclose(twiss_file);
    }

    const double I_1=R56;
    const double I_2=synch_rad_integral_2(line, args.periodicity);
    const double I_3=synch_rad_integral_3(line, args.periodicity);
    double I_4=0.0, I_5=0.0;
    for (size_t i=0; i<arrlenu(line); i++) {
      if (line[i].type == ELETYPE_SBEND) {
        double angle = line[i].as.sbend.angle;
        double L = line[i].as.sbend.length;
        double K1 = line[i].as.sbend.K1;
        double h = angle / L;
#if 0
        double eta = (element_etas[i]+element_etas[i+1]) / 2;
#else
        double (*sinlike_func)(double);
        double (*coslike_func)(double);
        double omega_sqr = pow(h, 2) + K1;
        double omega = sqrt(fabs(omega_sqr));
        double sign = 0.0;
        if (omega_sqr > 0.0) {
          sinlike_func = sin;
          coslike_func = cos;
          sign = 1.0;
        } else {
          sinlike_func = sinh;
          coslike_func = cosh;
          sign = -1.0;
        }
        double eta = element_etas[i] * sinlike_func(omega*L) / (omega*L);
        eta += sign * element_etaps[i] * (1 - coslike_func(omega*L)) / (omega*omega*L);
        eta += (h/K1) * (L - sinlike_func(omega*L)/omega);
#endif

        I_4 += args.periodicity * (eta * h * L * (2*K1 + h*h));
        I_5 += args.periodicity * (L * pow(fabs(h), 3) * element_curlyH[i]);
      }
    }

    const double gamma_0 = args.E_0 * 1e9 / ELECTRON_MASS;

    double j_x = 1.0f - I_4/I_2;
    double T_0 = total_length / C;

    printf("\nSynchrotron radiation integrals:\n");
    printf("\tI_1 = %+0.3e  (%+0.3e for the line)\n", I_1, I_1 / args.periodicity);
    printf("\tI_2 = %+0.3e  (%+0.3e for the line)\n", I_2, I_2 / args.periodicity);
    printf("\tI_3 = %+0.3e  (%+0.3e for the line)\n", I_3, I_3 / args.periodicity);
    printf("\tI_4 = %+0.3e  (%+0.3e for the line)\n", I_4, I_4 / args.periodicity);
    printf("\tI_5 = %+0.3e  (%+0.3e for the line)\n", I_5, I_5 / args.periodicity);
    printf("\n");
    printf("x fractional tune:    %f\n", acos(x_trace / 2) / (2*M_PI));
    printf("y fractional tune:    %f\n", acos(y_trace / 2) / (2*M_PI));
    printf("Energy loss per turn: %0.3f keV\n", e_loss_per_turn(I_2, gamma_0) / 1e3);
    printf("Momentum compaction:  %0.3e\n", R56 / total_length);
    printf("j_x:                  %+0.3e\n", j_x);
    printf("tau_x:                %0.3f ms\n", 
           1e3 * (2 * args.E_0*1e9 * T_0) / (j_x * e_loss_per_turn(I_2, gamma_0)));
    printf("Natural x emittance:  %0.3f pm.rad\n", 
           1e12 * natural_emittance_x(I_2, I_4, I_5, gamma_0));
    printf("Energy spread:        %0.3e\n", sqrt(energy_spread(I_2, I_3, I_4, gamma_0)));
  
    arrfree(element_etas);
    arrfree(element_etaps);
    arrfree(element_beta_xs);
    arrfree(element_beta_ys);
    arrfree(element_curlyH);
    arrfree(Ss);
  }

  printf("\n");

  arrfree(line);

  return 0;
}
