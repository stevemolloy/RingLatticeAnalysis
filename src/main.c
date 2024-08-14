#include <stdio.h>
#include <stdlib.h>

#include "lib.h"

int main(void)
{
  char *file_path = "./max4u_lattice.mad8";
  char *buffer = read_entire_file(file_path);

  printf("%s\n", buffer);

  free(buffer);

  return 0;
}
