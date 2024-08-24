#define _POSIX_C_SOURCE 200809L
#include <stdbool.h>

#include "lib.h"
#include "elements.h"

#define STB_DS_IMPLEMENTATION
#include "stb_ds.h"

#define ARENA_IMPLEMENTATION
#include "arena.h"

#define NOB_IMPLEMENTATION
#include "nob.h"

Element *element_list = NULL;
ElementLibrary *element_library = NULL;

int main(int argc, char **argv) {
  size_t periodicity = 1;
  size_t harmonic_number = 1;
  float E_0 = 0;
  char *file_path = NULL;

  nob_shift_args(&argc, &argv);
  while (argc > 0) {
    char *next_arg = nob_shift_args(&argc, &argv);
    if (strcmp(next_arg, "-p")==0) {
      periodicity = atoi(nob_shift_args(&argc, &argv));
    } else if (strcmp(next_arg, "-h")==0) {
      harmonic_number = atoi(nob_shift_args(&argc, &argv));
    } else if (strcmp(next_arg, "-E")==0) {
      E_0 = atoi(nob_shift_args(&argc, &argv));
    } else {
      file_path = next_arg;
    }
  }

  float gamma_0 = E_0 * 1e9 / ELECTRON_MASS;

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

  printf("\n");
  printf("Summary of the lattice defined in %s\n", file_path);
  printf("\n");
  printf("Periodicity: %zu\n", periodicity);
  printf("Harmonic number: %zu\n", harmonic_number);
  printf("Total length of the line: %f m\n", total_length);
  printf("Total length of %zu lines: %f m\n", periodicity, periodicity * total_length);
  printf("Total bending angle of the line: %0.3f degrees\n", radians_to_degrees(total_angle));
  printf("Total bending angle of %zu lines: %0.3f degrees\n", periodicity, periodicity * radians_to_degrees(total_angle));
  printf("RF frequency for a hamonic number of %zu: %0.6f MHz\n", harmonic_number, rf_freq/1e6);
  printf("\n");
  printf("Synchrotron radiation integrals:\n");
  printf("\tI_2 = %f\n", I2);
  printf("\tI_3 = %f\n", I3);
  printf("\n");
  printf("Energy loss per turn: %0.3f keV\n", e_loss_per_turn(I2, gamma_0) / 1e3);

  double mat1[] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
  double mat2[] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
  double res[] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
  matrix_multiply(mat1, mat2, res, 3, 3, 3, 3);
  for (size_t i=0; i<3; i++) {
    for (size_t j=0; j<3; j++) {
      printf("%f ", res[i*3 + j]);
    }
    printf("\n");
  }

  arrfree(line);
  arrfree(element_list);
  shfree(element_library);
  free(buffer);

  printf("\n");

  return 0;
}
