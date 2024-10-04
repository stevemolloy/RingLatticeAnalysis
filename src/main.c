#define _POSIX_C_SOURCE 200809L
#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "lib.h"
#include "elements.h"

#define M_PI acos(-1.0)

#define STB_DS_IMPLEMENTATION
#include "stb_ds.h"

#define NOB_IMPLEMENTATION
#include "nob.h"

Element *element_list = NULL;
ElementLibrary *element_library = NULL;

void usage(const char* program) {
  fprintf(stderr, "USAGE: %s [-p <periodicity>] [-h <harmonic number>] [-E <kinetic energy>] lattice_file\n", program);
}

int main(int argc, char **argv) {
  const char* program = nob_shift_args(&argc, &argv);

  int periodicity = 1;
  int harmonic_number = 1;
  double E_0 = 0;
  char *file_path = NULL;

  while (argc > 0) {
    char *next_arg = nob_shift_args(&argc, &argv);
    if (strcmp(next_arg, "-p")==0) {
      periodicity = atoi(nob_shift_args(&argc, &argv));
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

  double gamma_0 = E_0 * 1e9 / ELECTRON_MASS;

  char *buffer = read_entire_file(file_path);
  char *cursor = buffer;

  cursor = join_lines(cursor);

  cursor = populate_element_library(cursor);

  Element* line = NULL;
  create_line(cursor, &line);

  double total_length = calculate_line_length(line);
  double total_angle = calculate_line_angle(line);
  double rf_freq = C / (periodicity * total_length / harmonic_number);

  double I2 = synch_rad_integral_2(line, periodicity);
  double I3 = synch_rad_integral_3(line, periodicity);

  printf("Summary of the lattice defined in %s\n", file_path);
  printf("\n");
  printf("Periodicity: %d\n", periodicity);
  printf("Harmonic number: %d\n", harmonic_number);
  printf("Number of elements in the line: %td\n", arrlen(line));
  printf("Total length of the line: %f m\n", total_length);
  printf("Total length of %d lines: %f m\n", periodicity, periodicity * total_length);
  printf("Total bending angle of the line: %0.3f degrees\n", radians_to_degrees(total_angle));
  printf("Total bending angle of %d lines: %0.3f degrees\n", periodicity, periodicity * radians_to_degrees(total_angle));
  printf("RF frequency for a hamonic number of %d: %0.6f MHz\n", harmonic_number, rf_freq/1e6);
  printf("\n");
  printf("Synchrotron radiation integrals:\n");
  printf("\tI_2 = %f\n", I2);
  printf("\tI_3 = %f\n", I3);
  printf("\n");
  printf("Energy loss per turn: %0.3f keV\n", e_loss_per_turn(I2, gamma_0) / 1e3);

  // for (size_t i=0; i<arrlenu(line); i++) {
  //   if (line[i].type == ELETYPE_SBEND) {
  //     printf("\n");
  //     element_print(line[i]);
  //     rmatrix_print(line[i].R_matrix);
  //   }
  // }

  double line_matrix[BEAM_DOFS*BEAM_DOFS] = {
    1, 0, 0, 0, 0, 0,
    0, 1, 0, 0, 0, 0,
    0, 0, 1, 0, 0, 0,
    0, 0, 0, 1, 0, 0,
    0, 0, 0, 0, 1, 0,
    0, 0, 0, 0, 0, 1,
  };

  double total_matrix[BEAM_DOFS*BEAM_DOFS] = {
    1, 0, 0, 0, 0, 0,
    0, 1, 0, 0, 0, 0,
    0, 0, 1, 0, 0, 0,
    0, 0, 0, 1, 0, 0,
    0, 0, 0, 0, 1, 0,
    0, 0, 0, 0, 0, 1,
  };

  for (size_t i=0; i<arrlenu(line); i++) {
    double temp_result[BEAM_DOFS*BEAM_DOFS] = {0};
    matrix_multiply(line[i].R_matrix, line_matrix, temp_result, 6, 6, 6, 6);
    memcpy(line_matrix, temp_result, BEAM_DOFS*BEAM_DOFS*sizeof(double));
  }

  printf("\n");
  printf("Total matrix for the line is:\n");
  rmatrix_print(line_matrix);

  double x_trace = (line_matrix[0*BEAM_DOFS + 0] + line_matrix[1*BEAM_DOFS + 1]) / 2;
  double y_trace = (line_matrix[2*BEAM_DOFS + 2] + line_matrix[3*BEAM_DOFS + 3]) / 2;
  printf("\n");
  printf("x fractional tune = %f\n", acos(x_trace) / (2*M_PI));
  printf("y fractional tune = %f\n", acos(y_trace) / (2*M_PI));

  for (size_t i=0; i<(size_t)periodicity; i++) {
    double temp_result[BEAM_DOFS*BEAM_DOFS] = {0};
    matrix_multiply(line_matrix, total_matrix, temp_result, 6, 6, 6, 6);
    memcpy(total_matrix, temp_result, BEAM_DOFS*BEAM_DOFS*sizeof(double));
  }

  printf("\n");
  printf("Total matrix:\n");
  rmatrix_print(total_matrix);

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

  arrfree(line);
  arrfree(element_list);
  shfree(element_library);
  free(buffer);

  return 0;
}
