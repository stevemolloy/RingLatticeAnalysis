#define _POSIX_C_SOURCE 200809L
#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "lib.h"
#include "elements.h"

#define M_PI (acos(-1.0))
#define ANGLE_EPSILON (0.1*M_PI/180)

#include "stb_ds.h"

#define NOB_IMPLEMENTATION
#include "nob.h"

Element *element_list = NULL;
ElementLibrary *element_library = NULL;

double line_matrix[BEAM_DOFS*BEAM_DOFS] = {0};
double total_matrix[BEAM_DOFS*BEAM_DOFS] = {0};

void usage(const char* program) {
  fprintf(stderr, "USAGE: %s [-p <periodicity>] [-h <harmonic number>] [-E <kinetic energy>] lattice_file\n", program);
}

int main(int argc, char **argv) {
  const char* program = nob_shift_args(&argc, &argv);

  size_t periodicity = 1;
  int harmonic_number = 1;
  double E_0 = 0;
  char *file_path = NULL;

  while (argc > 0) {
    char *next_arg = nob_shift_args(&argc, &argv);
    if (strcmp(next_arg, "-p")==0) {
      periodicity = (size_t)atoi(nob_shift_args(&argc, &argv));
    } else if (strcmp(next_arg, "-h")==0) {
      harmonic_number = atoi(nob_shift_args(&argc, &argv));
    } else if (strcmp(next_arg, "-E")==0) {
      E_0 = strtod(nob_shift_args(&argc, &argv), NULL);
    } else {
      file_path = next_arg;
    }
  }

  if (file_path == NULL) {
    usage(program);
    return 1;
  }

  char *buffer = read_entire_file(file_path);
  char *cursor = buffer;

  cursor = join_lines(cursor);

  cursor = populate_element_library(cursor);

  Element* line = NULL;
  create_line(cursor, &line);

  printf("\nSummary of the lattice defined in %s\n", file_path);
  printf("\n");

  double line_length = calculate_line_length(line);
  double total_length = line_length * periodicity;

  double line_angle = calculate_line_angle(line);
  double total_angle = line_angle * periodicity;

  bool closed_system = !(fabs(total_angle - 2*M_PI) > ANGLE_EPSILON);
  if (!closed_system) {
    printf("Total bending angle (%0.3f degrees) is not within %0.1e degrees of a circle.\n", 
           radians_to_degrees(total_angle), radians_to_degrees(ANGLE_EPSILON));
    printf("System does not close, so not calculating ring parameters.\n");
    printf("\n");
  }

  get_line_matrix(line_matrix, line);
  apply_matrix_n_times(total_matrix, line_matrix, periodicity);

  printf("Periodicity: %zu\n", periodicity);
  printf("Harmonic number: %d\n", harmonic_number);
  printf("Number of elements in the line: %td\n", arrlen(line));
  printf("Total length of the line: %f m\n", line_length);
  if (periodicity != 1) {
    printf("Total length of %zu lines: %f m\n", periodicity, total_length);
  }
  printf("Total bending angle of the line: %0.3f degrees\n", radians_to_degrees(line_angle));
  if (periodicity != 1) {
    printf("Total bending angle of %zu lines: %0.3f degrees\n", 
           periodicity, radians_to_degrees(total_angle));
  }

  if (closed_system) {
    double rf_freq = C / (total_length / harmonic_number);

    double gamma_0 = E_0 * 1e9 / ELECTRON_MASS;
    double I2 = synch_rad_integral_2(line, periodicity);
    double I3 = synch_rad_integral_3(line, periodicity);

    printf("RF frequency for a hamonic number of %d: %0.6f MHz\n", harmonic_number, rf_freq/1e6);
    printf("\n");
    printf("Synchrotron radiation integrals:\n");
    printf("\tI_2 = %f\n", I2);
    printf("\tI_3 = %f\n", I3);
    printf("\n");
    printf("Energy loss per turn: %0.3f keV\n", e_loss_per_turn(I2, gamma_0) / 1e3);
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

  if (periodicity != 1) {
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

    double R11 = line_matrix[0*BEAM_DOFS + 0];
    double R12 = line_matrix[0*BEAM_DOFS + 1];
    double R21 = line_matrix[1*BEAM_DOFS + 0];
    double R22 = line_matrix[1*BEAM_DOFS + 1];
    double R16 = line_matrix[0*BEAM_DOFS + 5];
    double R26 = line_matrix[1*BEAM_DOFS + 5];

    double eta_x  = ((1 - R22)*R16 -    R12   *R26) / (1 - R11 - R22 + R11*R22 - R12*R21);
    double eta_px = (  -R21   *R16 + (1 - R11)*R26) / (1 - R11 - R22 + R11*R22 - R12*R21);

    printf("eta_x  = %f\n", eta_x);
    printf("eta_px = %f\n", eta_px);
  }

  printf("\n");

  arrfree(line);
  arrfree(element_list);
  shfree(element_library);
  free(buffer);

  return 0;
}
