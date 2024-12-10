#include <float.h>
#include <ctype.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "elements.h"
#include "lib.h"

#include "stb_ds.h"

bool compare_files(const char *testname, const char *filename1, const char *filename2);

bool test_matmul(void);
bool test_sbend(void);
bool test_full_lat_all_mats(void);
bool compare_with_matlab(void);
bool test_threebythree(void);
bool test_twiss_propagation(void);
bool test_synchrad_integrals(void);
bool test_populate_element_library(void);

typedef bool (*TestFunction)(void);

TestFunction test_functions[] = {
  test_matmul,
  test_sbend,
  /*test_full_lat_all_mats, compare_with_matlab*/
  test_threebythree,
  test_twiss_propagation,
  test_synchrad_integrals,
  test_populate_element_library,
};

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

bool test_threebythree(void) {
  const char *test_name = "MULTIPLY MATRIX BY VECTOR TEST";
  double threebythree[9] = {
    1, 2, 3,
    4, 5, 6,
    7, 8, 9,
  };
  double vec3[3] = {2, 2, 2};
  double result[3] = {0};
  double expected[3] = {12, 30, 48};

  matrix_multiply(threebythree, vec3, result, 3, 3, 3, 1);
  for (size_t i=0; i<3; i++) {
    if (result[i] != expected[i]) {
      printf("%s FAILED\n", test_name);
      return false;
    }
  }

  printf("%s PASSED\n", test_name);
  return true;
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
    -21.8411,  10.0408, -25.7030,  -1.1026, -29.2369,  49.5338,
    -49.3218,  28.1680, -15.7278,  30.6067, -11.7935,  20.7609,
     -0.4128, -38.8467,   4.5440, -12.2152,  16.0280, -41.9433,
     48.8482,   7.9329, -43.2427,   1.7978,  25.8373, -45.6692,
     23.7941,  37.0371,  -8.9552, -40.5402, -32.6931,  -0.8844,
    -18.9280,  18.9776, -26.2489,  40.9095,   1.7380,  -5.3404,
  };
  double matB[36] = {
     -1.3201,  40.7355, -21.0775, -26.0287,  -5.6128, -27.7119,
    -33.4109, -40.5693, -42.6913,  45.9483,  49.5819, -10.4391,
    -13.9343, -31.8674, -30.5392, -19.4538,  -6.3413, -27.5475,
     38.0722,  44.6587,  -8.2515, -34.5085, -19.5559, -22.9975,
     24.4352, -39.9151, -20.7073,   5.5508, -25.3490,  -8.1559,
     -8.3229, -11.1962,  20.2139,  29.0544,  46.0825,  49.7736,
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

  double matC[36] = {
    0.3171, 0.7952, 0.7547, 0.4984, 0.2551, 0.1386,
    0.9502, 0.1869, 0.2760, 0.9597, 0.5060, 0.1493,
    0.0344, 0.4898, 0.6797, 0.3404, 0.6991, 0.2575,
    0.4387, 0.4456, 0.6551, 0.5853, 0.8909, 0.8407,
    0.3816, 0.6463, 0.1626, 0.2238, 0.9593, 0.2543,
    0.7655, 0.7094, 0.1190, 0.7513, 0.5472, 0.8143,
  };
  double matD[36] = {
    0.2435, 0.4733, 0.2858, 0.0540, 0.4694, 0.5285,
    0.9293, 0.3517, 0.7572, 0.5308, 0.0119, 0.1656,
    0.3500, 0.8308, 0.7537, 0.7792, 0.3371, 0.6020,
    0.1966, 0.5853, 0.3804, 0.9340, 0.1622, 0.2630,
    0.2511, 0.5497, 0.5678, 0.1299, 0.7943, 0.6541,
    0.6160, 0.9172, 0.0759, 0.5688, 0.3112, 0.6892,
  };

  if (!matrix_multiply(matC, matD, result, 6, 6, 6, 6)) {
    fprintf(stderr, "Matrix multiplication didn't work");
    return false;
  };
  fprintf(result_file, "matC = \n");
  rmatrix_print(result_file, matC);
  fprintf(result_file, "\n");

  fprintf(result_file, "matD = \n");
  rmatrix_print(result_file, matD);
  fprintf(result_file, "\n");

  fprintf(result_file, "matC * matD = \n");
  rmatrix_print(result_file, result);
  fprintf(result_file, "\n");

  if (!matrix_multiply(matD, matC, result, 6, 6, 6, 6)) {
    fprintf(stderr, "Matrix multiplication didn't work");
    return false;
  };
  fprintf(result_file, "matD * matC = \n");
  rmatrix_print(result_file, result);
  fprintf(result_file, "\n");

  fclose(result_file);

  return compare_files(test_name, expected_filename, result_filename);
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

  arrfree(line);
  
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
  
  arrfree(line);

  return comparison_result;
}

bool test_synchrad_integrals(void) {
  const char *test_name = "SYNCH RAD INTEGRALS TEST";
  const char *expected_filename = "./tests/synchradintegrals_expected.txt";
  const char *result_filename =   "./tests/synchradintegrals_result.txt";

  char *file_path = "./lattices/max4_r3_lattice.mad8";

  Element *line = {0};
  generate_lattice(file_path, &line);

  const double periodicity = 20;

  double line_matrix[BEAM_DOFS*BEAM_DOFS] = {0};
  double total_matrix[BEAM_DOFS*BEAM_DOFS] = {0};

  get_line_matrix(line_matrix, line);
  apply_matrix_n_times(total_matrix, line_matrix, periodicity);

  LinOptsParams lin_opt_params = {
    .Ss = NULL,
    .element_beta_xs = NULL,
    .element_beta_ys = NULL,
    .element_etas = NULL,
    .element_etaps = NULL,
    .element_curlyH = NULL,
  };
  double I[5] = {0};
  propagate_linear_optics(line, line_matrix, &lin_opt_params, I);

  FILE *synradint_file = fopen(result_filename, "w");
  for (size_t i=0; i<5; i++) {
    fprintf(synradint_file, "I[%zu] = %e\n", i+1, I[i]);
  }
  fclose(synradint_file);

  return compare_files(test_name, expected_filename, result_filename);
}

bool test_populate_element_library(void) {
  const char *test_name = "ELEMENT PARSING TEST";
  const char *expected_filename = "./tests/eleparse_expected.txt";
  const char *result_filename =   "./tests/eleparse_result.txt";

  char *file_path = "./lattices/max4_r3_lattice.mad8";

  Element *line = {0};
  generate_lattice(file_path, &line);

  FILE *element_file = fopen(result_filename, "w");
  for (size_t i=0; i<arrlenu(line); i++) {
    element_print(element_file, line[i]);
  }
  fclose(element_file);

  return compare_files(test_name, expected_filename, result_filename);
}

bool test_twiss_propagation(void) {
  const char *test_name = "TWISS PROP TEST";
  const char *expected_filename = "./tests/twissprop_expected.txt";
  const char *result_filename =   "./tests/twissprop_result.txt";

  char *file_path = "./lattices/m4U_f02020101_lattice.mad8";

  Element *line = {0};
  generate_lattice(file_path, &line);

  const double periodicity = 20;

  double line_matrix[BEAM_DOFS*BEAM_DOFS] = {0};
  double total_matrix[BEAM_DOFS*BEAM_DOFS] = {0};

  LinOptsParams lin_opt_params = { 0 };
  double I[5] = {0};
  propagate_linear_optics(line, total_matrix, &lin_opt_params, I);

  apply_matrix_n_times(total_matrix, line_matrix, periodicity);

  FILE *twiss_file = fopen(result_filename, "w");
  fprintf(twiss_file, "S / m, beta_x / m, beta_y / m, eta_x / m\n");
  for (size_t i=0; i<arrlenu(line); i++) {
    fprintf(twiss_file, "%0.6e, %0.6e, %0.6e, %0.6e\n",
				lin_opt_params.Ss[i],
				lin_opt_params.element_beta_xs[i],
				lin_opt_params.element_beta_ys[i],
				lin_opt_params.element_etas[i]);
  }
  fclose(twiss_file);

  return compare_files(test_name, expected_filename, result_filename);
}

bool compare_arrays(double *correct_mat, double *compare_mat, size_t n) {
  bool result = true;
  for (size_t i=0; i<n; i++) {
    double absolute_diff = fabs(compare_mat[i] - correct_mat[i]);
    double relative_diff = correct_mat[i]==0.0 ? DBL_MAX : absolute_diff / correct_mat[i];
    if ((relative_diff > 1e-6) && (absolute_diff > 1e-9)) {
      printf("    fabs(%e - %e)\n", correct_mat[i], compare_mat[i]);
      printf("    abs diff = %e :: rel dif = %e\n", absolute_diff, relative_diff);
      result = false;
    }
  }
  return result;
}

bool compare_with_matlab(void) {
  const char *test_name = "MATLAB_COMPARISON TEST";

  const char *filename = "./lattices/m4U_240521_b03_03_07_06.mad8";

  Element *line = {0};
  generate_lattice(filename, &line);

  const char *matlab_output_filename = "./tests/matlab_output.txt";
  char *matlab_output_buffer = read_entire_file(matlab_output_filename);
  char *cursor = matlab_output_buffer;

  for (size_t i=0; i<arrlenu(line); i++) {
    if (strncmp(cursor, "Element ", 8) != 0) {
      printf("%s FAILED\n", test_name);
      printf("    Error while parsing element %zu\n", i+1);
      printf("    Expected the text \"Element\", but got something else. File not formatted correctly.\n");
      return false;
    }
    cursor += 8;
    size_t ele_num = strtol(cursor, &cursor, 10);
    if (ele_num != i+1) {
      printf("%s FAILED\n", test_name);
      printf("    Error while parsing element %zu\n", i+1);
      printf("    Got an unexpected element number from the matlab output file\n");
      return false;
    }
    while ((*cursor==':') || isspace(*cursor)) {
      cursor++;
    }
    if ((strncmp(cursor, "Drift", 5) == 0)) {
      cursor += 5;
      if (line[i].type != ELETYPE_DRIFT) {
        printf("%s FAILED\n", test_name);
        printf("    Error while parsing element %zu\n", i+1);
        printf("    Expected a DRIFT, but got something else.\n");
        return false;
      }
    } else if ((strncmp(cursor, "Marker", 6) == 0)) {
      cursor += 6;
      if (line[i].type != ELETYPE_DRIFT) {
        printf("%s FAILED\n", test_name);
        printf("    Error while parsing element %zu\n", i+1);
        printf("    Expected a DRIFT, but got something else.\n");
        return false;
      }
    } else if ((strncmp(cursor, "Monitor", 7) == 0)) {
      cursor += 7;
      if (line[i].type != ELETYPE_DRIFT) {
        printf("%s FAILED\n", test_name);
        printf("    Error while parsing element %zu\n", i+1);
        printf("    Expected a DRIFT, but got something else.\n");
        return false;
      }
    } else if ((strncmp(cursor, "Corrector", 9) == 0)) {
      cursor += 9;
      if (line[i].type != ELETYPE_DRIFT) {
        printf("%s FAILED\n", test_name);
        printf("    Error while parsing element %zu\n", i+1);
        printf("    Expected a DRIFT, but got something else.\n");
        return false;
      }
    } else if ((strncmp(cursor, "Multipole", 9) == 0)) {
      cursor += 9;
      if ((line[i].type != ELETYPE_SBEND) && 
        (line[i].type != ELETYPE_QUAD) && 
        (line[i].type != ELETYPE_OCTUPOLE) && 
        (line[i].type != ELETYPE_SEXTUPOLE) && 
        (line[i].type != ELETYPE_MULTIPOLE)) {
        printf("%s FAILED\n", test_name);
        printf("    Error while parsing element %zu\n", i+1);
        printf("    Expected something that matches a multipole, but got something else.\n");
        return false;
      }
    } else if ((strncmp(cursor, "Bend", 4) == 0)) {
      cursor += 4;
      if ((line[i].type != ELETYPE_SBEND)) {
        printf("%s FAILED\n", test_name);
        printf("    Error while parsing element %zu\n", i+1);
        printf("    Expected something that matches a multipole, but got something else.\n");
        return false;
      }
    } else if ((strncmp(cursor, "RFCavity", 8) == 0)) {
      cursor += 8;
      if ((line[i].type != ELETYPE_CAVITY)) {
        printf("%s FAILED\n", test_name);
        printf("    Error while parsing element %zu\n", i+1);
        printf("    Expected something that matches a multipole, but got something else.\n");
        return false;
      }
    } else {
      printf("%s FAILED\n", test_name);
      printf("    Error while parsing element %zu\n", i+1);
      printf("    Got an unknown element type.\n");
      return false;
    }

    double mat_from_matlab[BEAM_DOFS*BEAM_DOFS] = {0};
    for (size_t j=0; j<BEAM_DOFS*BEAM_DOFS; j++) {
      size_t y_ind = j / BEAM_DOFS;
      size_t x_ind = j % BEAM_DOFS;
      double sign = 1.0;
      if (y_ind == 4) {
        y_ind = 5;
        sign = -1.0;
      } else if (y_ind == 5) {
        y_ind = 4;
        sign = -1.0;
      }
      if (x_ind == 4) {
        x_ind = 5;
        sign = -1.0;
      } else if (x_ind == 5) {
        x_ind = 4;
        sign = -1.0;
      }
      mat_from_matlab[y_ind*BEAM_DOFS + x_ind] = sign * strtod(cursor, &cursor);
    }
    if (!compare_arrays(mat_from_matlab, line[i].R_matrix, BEAM_DOFS*BEAM_DOFS)) {
      printf("%s FAILED\n", test_name);
      printf("    Error while parsing element %zu\n", i+1);
      printf("    Matrices not equal.\n");
      return false;
    }

    while (isspace(*cursor)) {
      cursor++;
    }
  }

  free(matlab_output_buffer);
  arrfree(line);

  printf("%s PASSED\n", test_name);
  return true;
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

