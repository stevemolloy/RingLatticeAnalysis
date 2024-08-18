#include <stdbool.h>

#include "lib.h"
#include "elements.h"

#define STB_DS_IMPLEMENTATION
#include "stb_ds.h"

#define ARENA_IMPLEMENTATION
#include "arena.h"

Element *element_list = NULL;
ElementLibrary *element_library = NULL;

int main(void) {
  char *file_path = "./lattices/max4u_lattice.mad8";
  char *buffer = read_entire_file(file_path);
  char *cursor = buffer;

  cursor = join_lines(cursor);

  cursor = populate_element_library(cursor);

  Element* line = NULL;
  create_line(cursor, &line);

  float total_length = calculate_line_length(line);

  printf("\n");
  printf("Total length of this line is %f m\n", total_length);
  printf("Total length of 20 lines is %f m\n", 20 * total_length);

  int harmonic_number = 176;
  float rf_freq = C / (20 * total_length / harmonic_number);
  printf("RF frequency for a hamonic number of %i is %0.3f MHz\n", harmonic_number, rf_freq/1e6);

  arrfree(line);
  arrfree(element_list);
  shfree(element_library);
  free(buffer);

  printf("\n");

  return 0;
}
