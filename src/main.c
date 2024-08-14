#include <stdio.h>
#include <stdlib.h>

#include "lib.h"

typedef struct {
  float length;
} Drift;

int main(void) {
  char *file_path = "./lattices/max4u_lattice.mad8";
  char *buffer = read_entire_file(file_path);

  char *cursor = buffer;

  while (*cursor != '\0') {
    if (*cursor == '!') {
      advance_to_next_line(&cursor);
    } else {
      printf("%c", *cursor);
      cursor++;
    }
  }

  free(buffer);

  return 0;
}
