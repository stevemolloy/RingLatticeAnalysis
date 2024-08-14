#include <stdio.h>

#include "lib.h"

int main(void)
{
  char *file_path = "./max_4u_f_0_20240807.m";
  char *buffer;
  char **lines;
  size_t num_lines = read_entire_file_to_lines(file_path, &buffer, &lines);

  printf("I read %zu lines\n", num_lines);

  return 0;
}
