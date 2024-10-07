#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "elements.h"
#include "lib.h"

Element *element_list = NULL;
ElementLibrary *element_library = NULL;

bool test_matmul(void);

int main(void) {
  bool result = true;
  size_t test_counter = 0;
  size_t passed_counter = 0;

  test_counter++;
  if (test_matmul()) passed_counter++;
  else result = false;

  printf("\n");
  printf("========================\n");
  printf("Tests: %zu (Passed: %zu)\n", test_counter, passed_counter);
  printf("========================\n");

  if (!result) return 1;
  return 0;
}

bool test_matmul(void) {
  const char *test_name = "MATRIX MULTIPLICATION TEST";
  const char *expected_filename = "./tests/matmul_expected.txt";
  const char *result_filename =   "./tests/matmul_result.txt";

  FILE *result_file = fopen(result_filename, "w");
  double matA[36] = {
    0.3171, 0.7952, 0.7547, 0.4984, 0.2551, 0.1386,
    0.9502, 0.1869, 0.2760, 0.9597, 0.5060, 0.1493,
    0.0344, 0.4898, 0.6797, 0.3404, 0.6991, 0.2575,
    0.4387, 0.4456, 0.6551, 0.5853, 0.8909, 0.8407,
    0.3816, 0.6463, 0.1626, 0.2238, 0.9593, 0.2543,
    0.7655, 0.7094, 0.1190, 0.7513, 0.5472, 0.8143,
  };
  double matB[36] = {
    0.2435, 0.4733, 0.2858, 0.0540, 0.4694, 0.5285,
    0.9293, 0.3517, 0.7572, 0.5308, 0.0119, 0.1656,
    0.3500, 0.8308, 0.7537, 0.7792, 0.3371, 0.6020,
    0.1966, 0.5853, 0.3804, 0.9340, 0.1622, 0.2630,
    0.2511, 0.5497, 0.5678, 0.1299, 0.7943, 0.6541,
    0.6160, 0.9172, 0.0759, 0.5688, 0.3112, 0.6892,
  };
  double result[36] = {0};
  if (!matrix_multiply(matA, matB, result, 6, 6, 6, 6)) {
    fprintf(stderr, "Matrix multiplication didn't work");
    return 1;
  };
  fprintf(result_file, "matA = \n");
  rmatrix_print(result_file, matA);
  fprintf(result_file, "\n");

  fprintf(result_file, "matB = \n");
  rmatrix_print(result_file, matB);
  fprintf(result_file, "\n");

  fprintf(result_file, "matA * matB = \n");
  rmatrix_print(result_file, result);
  fprintf(result_file, "\n");

  if (!matrix_multiply(matB, matA, result, 6, 6, 6, 6)) {
    fprintf(stderr, "Matrix multiplication didn't work");
    return 1;
  };
  fprintf(result_file, "matB * matA = \n");
  rmatrix_print(result_file, result);
  fprintf(result_file, "\n");

  fclose(result_file);

  char *expected_buff = read_entire_file(expected_filename);
  char *result_buff = read_entire_file(result_filename);
  
  bool comparison_result = true;
  if (strcmp(expected_buff, result_buff) != 0) {
    printf("%s FAILED\n", test_name);
    printf("    Expected (%s) does not match actual (%s)\n", expected_filename, result_filename);
    comparison_result = false;
  } else {
    printf("%s PASSED\n", test_name);
  }

  free(expected_buff);
  free(result_buff);

  return comparison_result;
}

