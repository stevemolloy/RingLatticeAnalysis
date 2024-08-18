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

  float total_length = calculate_line_length(line);
  float rf_freq = C / (periodicity * total_length / harmonic_number);

  printf("\n");
  printf("Total length of this line is %f m\n", total_length);
  printf("Total length of %zu lines is %f m\n", periodicity, periodicity * total_length);
  printf("RF frequency for a hamonic number of %zu is %0.3f MHz\n", harmonic_number, rf_freq/1e6);
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
