#include <stdbool.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "lib.h"
#include "elements.h"

#define M_PI (acos(-1.0))
#define ANGLE_EPSILON (0.1*M_PI/180)

#include "stb_ds.h"

void usage(const char* program) {
  fprintf(stderr, "USAGE: %s [-p <periodicity>] [-h <harmonic number>] [-E <kinetic energy>] lattice_file\n", program);
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
  double total_matrix[BEAM_DOFS*BEAM_DOFS] = {0};

  // for (size_t i=0; i<arrlenu(line); i++) {
  //   printf("%s :: ", line[i].name);
  //   switch (line[i].type) {
  //     case ELETYPE_DRIFT:     printf("DRIFT\n"); break;
  //     case ELETYPE_QUAD:      printf("QUAD\n"); break;
  //     case ELETYPE_CAVITY:    printf("CAVITY\n"); break;
  //     case ELETYPE_SEXTUPOLE: printf("SEXTUPOLE\n"); break;
  //     case ELETYPE_OCTUPOLE:  printf("OCTUPOLE\n"); break;
  //     case ELETYPE_MULTIPOLE: printf("MULTIPOLE\n"); break;
  //     case ELETYPE_SBEND: {
  //       printf("SBEND :: h = %0.16e\n", line[i].as.sbend.angle / line[i].as.sbend.length);
  //     }
  //       break;
  //   }
  //   rmatrix_print(stdout, line[i].R_matrix);
  //   printf("\n");
  // }
  get_line_matrix(line_matrix, line);

  printf("\nSummary of the lattice defined in %s\n", args.file_path);
  printf("\n");

  bool closed_system = !(fabs(total_angle - 2*M_PI) > ANGLE_EPSILON);
  if (!closed_system) {
    printf("Total bending angle (%0.3f degrees) is not within %0.1e degrees of a circle.\n", 
           radians_to_degrees(total_angle), radians_to_degrees(ANGLE_EPSILON));
    printf("System does not close, so not calculating ring parameters.\n");
    printf("\n");
  }

  printf("Periodicity: %zu\n", args.periodicity);
  printf("Harmonic number: %d\n", args.harmonic_number);
  printf("Number of elements in the line: %td\n", arrlen(line));
  printf("Total length of the line: %f m\n", line_length);
  if (args.periodicity != 1) {
    printf("Total length of %zu lines: %f m\n", args.periodicity, total_length);
  }
  printf("Total bending angle of the line: %0.3f degrees\n", radians_to_degrees(line_angle));
  if (args.periodicity != 1) {
    printf("Total bending angle of %zu lines: %0.3f degrees\n", 
           args.periodicity, radians_to_degrees(total_angle));
  }

  double I_1, I_2, I_3, I_4;//, I_5;
  if (closed_system) {
    double rf_freq = C / (total_length / args.harmonic_number);

    double gamma_0 = args.E_0 * 1e9 / ELECTRON_MASS;
    I_2 = synch_rad_integral_2(line, args.periodicity);
    I_3 = synch_rad_integral_3(line, args.periodicity);

    printf("RF frequency for a hamonic number of %d: %0.6f MHz\n", args.harmonic_number, rf_freq/1e6);
    printf("\n");
    printf("Energy loss per turn: %0.3f keV\n", e_loss_per_turn(I_2, gamma_0) / 1e3);
  }

  printf("\n");
  printf("Total matrix, R, for the line is:\n");
  rmatrix_print(stdout, line_matrix);

  double x_trace, y_trace;
  if (closed_system) {
    x_trace = (line_matrix[0*BEAM_DOFS + 0] + line_matrix[1*BEAM_DOFS + 1]) / 2;
    y_trace = (line_matrix[2*BEAM_DOFS + 2] + line_matrix[3*BEAM_DOFS + 3]) / 2;
    printf("\n");
    printf("x fractional tune = %f\n", acos(x_trace) / (2*M_PI));
    printf("y fractional tune = %f\n", acos(y_trace) / (2*M_PI));
  }

  if (args.periodicity != 1) {
    apply_matrix_n_times(total_matrix, line_matrix, args.periodicity);
    printf("\n");
    printf("Total matrix, R, for the full system:\n");
    rmatrix_print(stdout, total_matrix);
  }

  if (closed_system) {
    x_trace = (total_matrix[0*BEAM_DOFS + 0] + total_matrix[1*BEAM_DOFS + 1]) / 2;
    y_trace = (total_matrix[2*BEAM_DOFS + 2] + total_matrix[3*BEAM_DOFS + 3]) / 2;
    printf("\n");
    printf("x fractional tune = %f\n", acos(x_trace) / (2*M_PI));
    printf("y fractional tune = %f\n", acos(y_trace) / (2*M_PI));

    const double R11 = total_matrix[0*BEAM_DOFS + 0];
    const double R16 = total_matrix[0*BEAM_DOFS + 5];

    const double eta_x  = R16 / (1 - R11);

    printf("\n");
    printf("eta_x at beginning of the line  = %0.6e\n", eta_x);

    double *element_etas = {0};
    double *element_etaps = {0};
    double eta_vec[3] = {eta_x, 0.0f, 1.0f};
    arrput(element_etas, eta_vec[0]);
    arrput(element_etaps, eta_vec[1]);
    for (size_t i=0; i<arrlenu(line); i++) {
      double temp_vec[3] = {0};
      matrix_multiply(line[i].eta_prop_matrix, eta_vec, temp_vec, 3, 3, 3, 1);
      memcpy(eta_vec, temp_vec, 3*sizeof(double));
      arrput(element_etas, eta_vec[0]);
      arrput(element_etaps, eta_vec[1]);
    }

    I_1 = 0.0f;
    I_4 = 0.0f;
    for (size_t i=0; i<arrlenu(line); i++) {
      if (line[i].type == ELETYPE_SBEND) {
        double eta = (element_etas[i]+element_etas[i+1]) / 2;
        double h = line[i].as.sbend.angle / line[i].as.sbend.length;
        double K = line[i].as.sbend.K1;
        double L = line[i].as.sbend.length;

        I_1 += h * L * eta;
        I_4 += eta * h * L * (2*K + h*h);
      }
    }

    I_1 *= args.periodicity;
    I_4 *= args.periodicity;

    printf("\n");
    printf("Synchrotron radiation integrals:\n");
    printf("\tI_1 = %0.3e  (%0.3e for the line)\n", I_1, I_1 / args.periodicity);
    printf("\tI_2 = %0.3e  (%0.3e for the line)\n", I_2, I_2 / args.periodicity);
    printf("\tI_3 = %0.3e  (%0.3e for the line)\n", I_3, I_3 / args.periodicity);
    printf("\tI_4 = %0.3e  (%0.3e for the line (should be -2.049e-02))\n", I_4, I_4 / args.periodicity);
    printf("\n");
    printf("Momentum compaction = %0.3e\n", I_1 / total_length);
  
    arrfree(element_etas);
    arrfree(element_etaps);
  }

  printf("\n");

  arrfree(line);

  return 0;
}
