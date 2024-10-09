#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "elements.h"
#include "lib.h"

#include "stb_ds.h"

bool test_matmul(void);
bool test_sbend(void);
bool test_full_lat_all_mats(void);

typedef bool (*TestFunction)(void);

TestFunction test_functions[] = {test_matmul, test_sbend, test_full_lat_all_mats};

int main(void) {
  bool result = true;
  size_t test_counter = 0;
  size_t passed_counter = 0;
  size_t num_of_functions = (sizeof(test_functions)/sizeof(test_functions[0]));

  for (size_t i=0; i<num_of_functions; i++) {
    test_counter++;
    if (test_functions[i]()) passed_counter++;
    else result = false;
  }

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

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
  FILE *result_file;
  fopen_s(&result_file, result_filename, "w");
#else
  FILE *result_file = fopen(result_filename, "w");
#endif
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
    return false;
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
    return false;
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

bool compare_files(const char *testname, const char *filename1, const char *filename2) {
  bool comparison_result = true;
  char *expected_buff = read_entire_file(filename1);
  char *result_buff = read_entire_file(filename2);
  
  if (strcmp(expected_buff, result_buff) != 0) {
    comparison_result = false;
  } else {
    comparison_result = true;
  }

  if (!comparison_result) {
    printf("%s FAILED\n", testname);
    printf("    Expected (%s) does not match actual (%s)\n", filename1, filename2);
    comparison_result = false;
  } else {
    printf("%s PASSED\n", testname);
  }

  free(expected_buff);
  free(result_buff);

  return comparison_result;
}

bool test_sbend(void) {
  const char *test_name = "SBEND ELEMENT TEST";
  const char *expected_filename = "./tests/sbend_expected.txt";
  const char *result_filename =   "./tests/sbend_result.txt";

  const char *filename = "./lattices/whiskey.mad8";

  Element *line = {0};
  generate_lattice(filename, &line);
  double line_matrix[BEAM_DOFS*BEAM_DOFS] = {0};
  get_line_matrix(line_matrix, line);

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
  FILE *result_file;
  fopen_s(&result_file, result_filename, "w");
#else
  FILE *result_file = fopen(result_filename, "w");
#endif

  fprintf(result_file, "Total matrix, R, for the line is:\n");
  rmatrix_print(result_file, line_matrix);

  fclose(result_file);

  bool comparison_result = compare_files(test_name, expected_filename, result_filename);
  
  return comparison_result;
}

bool test_full_lat_all_mats(void) {
  const char *test_name = "FULL LAT ALL MATS TEST";
  const char *expected_filename = "./tests/fulllat_allmats_expected.txt";
  const char *result_filename =   "./tests/fulllat_allmats_result.txt";

  const char *filename = "./lattices/m4U_240521_b03_03_07_06.mad8";

  Element *line = {0};
  generate_lattice(filename, &line);

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
  FILE *result_file;
  fopen_s(&result_file, result_filename, "w");
#else
  FILE *result_file = fopen(result_filename, "w");
#endif

  for (size_t i=0; i<arrlenu(line); i++) {
    fprintf(result_file, "Element %03zu: ", i+1);
    switch (line[i].type) {
      case ELETYPE_DRIFT:     fprintf(result_file, "DRIFT\n"); break;
      case ELETYPE_QUAD:      fprintf(result_file, "QUAD\n"); break;
      case ELETYPE_SBEND:     fprintf(result_file, "SBEND\n"); break;
      case ELETYPE_CAVITY:    fprintf(result_file, "CAVITY\n"); break;
      case ELETYPE_SEXTUPOLE: fprintf(result_file, "SEXTUPOLE\n"); break;
      case ELETYPE_OCTUPOLE:  fprintf(result_file, "OCTUPOLE\n"); break;
      case ELETYPE_MULTIPOLE: fprintf(result_file, "MULTIPOLE\n"); break;
    }
    rmatrix_print(result_file, line[i].R_matrix);
    fprintf(result_file, "\n");
  }

  fclose(result_file);

  bool comparison_result = compare_files(test_name, expected_filename, result_filename);
  
  return comparison_result;
}

