#include <stdbool.h>

#include "lib.h"
#include "elements.h"

#define STB_DS_IMPLEMENTATION
#include "stb_ds.h"

#define ARENA_IMPLEMENTATION
#include "arena.h"

Element *element_list = NULL;
ElementLibrary *element_library = NULL;

int main(int argc, char **argv) {
  size_t periodicity = 1;
  size_t harmonic_number = 1;
  char *file_path = NULL;

  argc--;
  argv++;
  while (argc > 0) {
    if (strcmp(*argv, "-p")==0) {
      argc--;
      argv++;
      periodicity = atoi(*argv);
    } else if (strcmp(*argv, "-h")==0) {
      argc--;
      argv++;
      harmonic_number = atoi(*argv);
    } else {
      file_path = *argv;
    }
    argc--;
    argv++;
  }

  char *buffer = read_entire_file(file_path);
  char *cursor = buffer;

  cursor = join_lines(cursor);

  cursor = populate_element_library(cursor);

  Element* line = NULL;
  create_line(cursor, &line);

  double total_length = calculate_line_length(line);
  double total_angle = calculate_line_angle(line);
  double rf_freq = C / (periodicity * total_length / harmonic_number);

  printf("\n");
  printf("Summary of the lattice defined in %s\n", file_path);
  printf("Periodicity: %zu\n", periodicity);
  printf("Harmonic number: %zu\n", harmonic_number);
  printf("Total length of the line: %f m\n", total_length);
  printf("Total length of %zu lines: %f m\n", periodicity, periodicity * total_length);
  printf("Total bending angle of the line: %0.3f degrees\n", radians_to_degrees(total_angle));
  printf("Total bending angle of %zu lines: %0.3f degrees\n", periodicity, periodicity * radians_to_degrees(total_angle));
  printf("RF frequency for a hamonic number of %zu: %0.6f MHz\n", harmonic_number, rf_freq/1e6);
  printf("Synchrotron radiation integrals:\n");
  printf("\tI_2 = %f\n", synch_rad_integral_2(line, periodicity));
  printf("\tI_3 = %f\n", synch_rad_integral_3(line, periodicity));

  arrfree(line);
  arrfree(element_list);
  shfree(element_library);
  free(buffer);

  printf("\n");

  return 0;
}
