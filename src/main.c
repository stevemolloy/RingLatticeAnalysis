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

#include "sdm_lib.h"

sdm_arena_t mem_arena = {0};

// TODO: Factor out the twiss propagation stuff into diagonalised 2x2 instead of the 3x3 stuff
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
  if (args.error) {
    usage(args.programname);
    return 1;
  }

  sdm_arena_init(&mem_arena, SDM_ARENA_DEFAULT_CAP);

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

  if (args.E_0 == 0) 
    printf("WARNING: Kinetic energy not provided, so not calculating all parameters\n\n");
  printf("Periodicity: %zu\n", args.periodicity);
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

  printf("\nTotal matrix, R, for the line is:\n");
  rmatrix_print(stdout, line_matrix);

  double total_matrix[BEAM_DOFS*BEAM_DOFS] = {0};
  apply_matrix_n_times(total_matrix, line_matrix, args.periodicity);
  if (args.periodicity > 1) {
    printf("\nTotal matrix, R, for the full system:\n");
    rmatrix_print(stdout, total_matrix);
  }

  if (closed_system) {
    double x_trace = (total_matrix[0*BEAM_DOFS + 0] + total_matrix[1*BEAM_DOFS + 1]);
    double y_trace = (total_matrix[2*BEAM_DOFS + 2] + total_matrix[3*BEAM_DOFS + 3]);

    const double R56 = total_matrix[4*BEAM_DOFS + 5];

    LinOptsParams lin_opt_params = {
      .Ss = NULL,
      .element_beta_xs = NULL,
      .element_beta_ys = NULL,
      .element_etas = NULL,
      .element_etaps = NULL,
      .element_curlyH = NULL,
    };
    propagate_linear_optics(line, total_matrix, &lin_opt_params);

    if (args.save_twiss) {
      FILE *twiss_file = fopen(args.twiss_filename, "w");
      fprintf(twiss_file, "S / m, beta_x / m, beta_y / m, eta_x / m\n");
      for (size_t i=0; i<arrlenu(line); i++) {
        fprintf(twiss_file, "%0.6e, %0.6e, %0.6e, %0.6e\n",
		    		lin_opt_params.Ss[i],
		    		lin_opt_params.element_beta_xs[i],
		    		lin_opt_params.element_beta_ys[i],
		    		lin_opt_params.element_etas[i]);
      }
      fclose(twiss_file);
    }

    const double I_1=R56;
    const double I_2=synch_rad_integral_2(line, args.periodicity);
    const double I_3=synch_rad_integral_3(line, args.periodicity);
    double I_4=0.0, I_5=0.0;
    for (size_t i=0; i<arrlenu(line); i++) {
      // if (line[i].type == ELETYPE_SBEND) {
      //   element_print(line[i]);
      // }
      if (line[i].type == ELETYPE_SBEND) {
        double angle = line[i].as.sbend.angle;
        double L = line[i].as.sbend.length;
        double K1 = line[i].as.sbend.K1;
        double h = angle / L;

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
        double mean_eta = lin_opt_params.element_etas[i] * sinlike_func(omega*L) / (omega*L)
                  + sign * lin_opt_params.element_etaps[i] * (1 - coslike_func(omega*L)) / (omega*omega*L)
                  + sign* h * (omega*L - sinlike_func(omega*L))/(pow(omega,3)*L);

        I_4 += args.periodicity * (mean_eta * h * L * (2*K1 + h*h));
        I_5 += args.periodicity * (L * pow(fabs(h), 3) * lin_opt_params.element_curlyH[i]);
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
    if (gamma_0 > 0)
      printf("Energy loss per turn: %0.3f keV\n", e_loss_per_turn(I_2, gamma_0) / 1e3);
    printf("Momentum compaction:  %0.3e\n", R56 / total_length);
    printf("j_x:                  %+0.3e\n", j_x);
    if (gamma_0 > 0) {
      printf("tau_x:                %0.3f ms\n", 
           1e3 * (2 * args.E_0*1e9 * T_0) / (j_x * e_loss_per_turn(I_2, gamma_0)));
      printf("Natural x emittance:  %0.3f pm.rad\n", 
           1e12 * natural_emittance_x(I_2, I_4, I_5, gamma_0));
      printf("Energy spread:        %0.3e\n", sqrt(energy_spread(I_2, I_3, I_4, gamma_0)));
    }

    assert(arrlenu(line) == arrlenu(lin_opt_params.element_curlyH));

    FILE *curlyhfile = fopen("curlyH.csv", "w");
    for (size_t i=0; i<arrlenu(lin_opt_params.element_curlyH); i++) {
      fprintf(curlyhfile, "%e, %e\n", element_length(line[i]), lin_opt_params.element_curlyH[i]);
    }
    fclose(curlyhfile);

    arrfree(lin_opt_params.element_etas);
    arrfree(lin_opt_params.element_etaps);
    arrfree(lin_opt_params.element_beta_xs);
    arrfree(lin_opt_params.element_beta_ys);
    arrfree(lin_opt_params.element_curlyH);
    arrfree(lin_opt_params.Ss);
  }

  // for (size_t i=0; i<arrlenu(line); i++) {
  //   element_print(line[i]);
  //   rmatrix_print(stdout, line[i].R_matrix);
  //   printf("\n");
  // }

  printf("\n");

  arrfree(line);

  return 0;
}
