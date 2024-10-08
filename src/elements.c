#define _GNU_SOURCE
#include <math.h>
#include <assert.h>
#include <string.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>

#define STB_DS_IMPLEMENTATION
#include "stb_ds.h"

#include "elements.h"
#include "lib.h"

#define EPSILON 1e-9

extern Element* element_list;
extern ElementLibrary *element_library;

static void calc_sbend_matrix(Element *element);
static void calc_quad_matrix(Element *element);
static void make_r_matrix(Element *element);
static void getMinor(double* matrix, double* minor, int n, int r, int c);

double synch_rad_integral_1(Element *line, int periodicity) {
  (void)line;
  (void)periodicity;
  return 0;
}

double synch_rad_integral_2(Element *line, int periodicity) {
  double I_2 = 0;
  for (size_t i=0; i<arrlenu(line); i++) {
    double rho = bending_radius_of_element(line[i]);
    I_2 += element_length(line[i]) / pow(rho, 2);
  }

  return I_2 * (double)periodicity;
}

double synch_rad_integral_3(Element *line, int periodicity) {
  double I_3 = 0;
  for (size_t i=0; i<arrlenu(line); i++) {
    double rho_abs = fabs(bending_radius_of_element(line[i]));
    I_3 += element_length(line[i]) / pow(rho_abs, 3);
  }

  return I_3 * (double)periodicity;
}

double element_length(Element element) {
  switch (element.type) {
    case ELETYPE_SBEND:
      return element.as.sbend.length;
    case ELETYPE_DRIFT:
      return element.as.drift.length;
    case ELETYPE_QUAD:
      return element.as.quad.length;
    case ELETYPE_SEXTUPOLE:
      return element.as.sextupole.length;
    case ELETYPE_OCTUPOLE:
      return element.as.octupole.length;
    case ELETYPE_MULTIPOLE:
      return element.as.multipole.length;
    case ELETYPE_CAVITY:
      return element.as.cavity.length;
  }
}

double bending_radius_of_element(Element element) {
  double rho = DBL_MAX;
  switch (element.type) {
    case ELETYPE_SBEND:
      if (element.as.sbend.angle != 0) {
        rho = element.as.sbend.length / element.as.sbend.angle;
      }
      break;
    case ELETYPE_DRIFT:
    case ELETYPE_QUAD:
    case ELETYPE_SEXTUPOLE:
    case ELETYPE_OCTUPOLE:
    case ELETYPE_MULTIPOLE:
    case ELETYPE_CAVITY:
      break;
  }
  return rho;
}

static void make_r_matrix(Element *element) {
  double ele_length = element_length(*element);

  memset(element->R_matrix, 0, sizeof(element->R_matrix));
  for (size_t i=0; i<BEAM_DOFS; i++) {
    element->R_matrix[i*BEAM_DOFS + i] = 1.0;
  }

  switch (element->type) {
    case (ELETYPE_DRIFT):
    case (ELETYPE_CAVITY):
    case (ELETYPE_MULTIPOLE):
    case (ELETYPE_SEXTUPOLE):
    case (ELETYPE_OCTUPOLE):
      element->R_matrix[0*BEAM_DOFS + 1] = ele_length;
      element->R_matrix[2*BEAM_DOFS + 3] = ele_length;
      break;
    case (ELETYPE_QUAD):
      calc_quad_matrix(element);
      break;
    case (ELETYPE_SBEND):
      calc_sbend_matrix(element);
      break;
  }
}

static void calc_quad_matrix(Element *element) {
  assert(element->type == ELETYPE_QUAD && "This func should only be called for quads.");
  double omega = sqrt(fabs(element->as.quad.K1));
  double L = element->as.quad.length;
  double omega_L = omega * L;

  if (element->as.quad.K1 == 0) {
    element->R_matrix[0*BEAM_DOFS + 1] = L;
    element->R_matrix[0*BEAM_DOFS + 1] = L;
    return;
  }

  if (element->as.quad.K1 > 0) {
    // Focusing 2x2
    element->R_matrix[0*BEAM_DOFS + 0] = cos(omega_L);
    element->R_matrix[0*BEAM_DOFS + 1] = sin(omega_L) / omega;
    element->R_matrix[1*BEAM_DOFS + 0] = sin(omega_L) * (-omega);
    element->R_matrix[1*BEAM_DOFS + 1] = cos(omega_L);
    // Defocusing 2x2
    element->R_matrix[2*BEAM_DOFS + 2] = cosh(omega_L);
    element->R_matrix[2*BEAM_DOFS + 3] = sinh(omega_L) / omega;
    element->R_matrix[3*BEAM_DOFS + 2] = sinh(omega_L) * omega;
    element->R_matrix[3*BEAM_DOFS + 3] = cosh(omega_L);
  } else if (element->as.quad.K1 < 0) {
    // Defocusing 2x2
    element->R_matrix[2*BEAM_DOFS + 2] = cos(omega_L);
    element->R_matrix[2*BEAM_DOFS + 3] = sin(omega_L) / omega;
    element->R_matrix[3*BEAM_DOFS + 2] = sin(omega_L) * (-omega);
    element->R_matrix[3*BEAM_DOFS + 3] = cos(omega_L);
    // Focusing 2x2
    element->R_matrix[0*BEAM_DOFS + 0] = cosh(omega_L);
    element->R_matrix[0*BEAM_DOFS + 1] = sinh(omega_L) / omega;
    element->R_matrix[1*BEAM_DOFS + 0] = sinh(omega_L) * omega;
    element->R_matrix[1*BEAM_DOFS + 1] = cosh(omega_L);
  }
}

static void calc_sbend_matrix(Element *element) {
  assert(element->type == ELETYPE_SBEND && "This func should only be called for sbends.");
  double K1 = element->as.sbend.K1;
  double L = element->as.sbend.length;
  double angle = element->as.sbend.angle;
  double rho = L / angle;
  double h = 1 / rho;

  double omega_x_sqr = pow(h, 2) + K1;
  double omega_x = sqrt(fabs(omega_x_sqr));
  double omega_x_L = omega_x * L;

  double omega_y_sqr = K1;
  double omega_y = sqrt(fabs(omega_y_sqr));
  double omega_y_L = omega_y * L;

  double (*sinlike_func)(double);
  double (*coslike_func)(double);
  double sign;

  if (omega_x == 0.0f) {
    element->R_matrix[0*BEAM_DOFS + 0] = 1;
    element->R_matrix[0*BEAM_DOFS + 1] = L;
    element->R_matrix[1*BEAM_DOFS + 0] = 0;
    element->R_matrix[1*BEAM_DOFS + 1] = 1;

    element->R_matrix[4*BEAM_DOFS + 4] = 1;
    element->R_matrix[4*BEAM_DOFS + 5] = 0;
    element->R_matrix[5*BEAM_DOFS + 4] = 0;
    element->R_matrix[5*BEAM_DOFS + 5] = 1;

    element->R_matrix[4*BEAM_DOFS + 0] = h * L;
    element->R_matrix[4*BEAM_DOFS + 1] = 0;
  } else {
    if (omega_x_sqr > 0) {
      sinlike_func = &sin;
      coslike_func = &cos;
      sign = -1;
    } else {
      sinlike_func = &sinh;
      coslike_func = &cosh;
      sign = 1;
    }
    assert((sign == 1) | (sign == -1) && "Sign value has been corrupted");
    element->R_matrix[0*BEAM_DOFS + 0] = coslike_func(omega_x_L);
    element->R_matrix[0*BEAM_DOFS + 1] = (1/omega_x) * sinlike_func(omega_x_L);
    element->R_matrix[1*BEAM_DOFS + 0] = sign * omega_x * sinlike_func(omega_x_L);
    element->R_matrix[1*BEAM_DOFS + 1] = coslike_func(omega_x_L);

    element->R_matrix[0*BEAM_DOFS + 5] = (h/pow(omega_x,2)) * (1 - coslike_func(omega_x_L));
    element->R_matrix[1*BEAM_DOFS + 5] = (h/omega_x) * sinlike_func(omega_x_L);

    element->R_matrix[4*BEAM_DOFS + 0] = -element->R_matrix[1*BEAM_DOFS + 5];
    element->R_matrix[4*BEAM_DOFS + 1] = -element->R_matrix[0*BEAM_DOFS + 5];

    element->R_matrix[4*BEAM_DOFS + 4] = 1;
    element->R_matrix[4*BEAM_DOFS + 5] = -pow(h,2)*(omega_x_L - sinlike_func(omega_x_L)) / pow(omega_x, 3);
    element->R_matrix[5*BEAM_DOFS + 4] = 0;
    element->R_matrix[5*BEAM_DOFS + 5] = 1;
  }

  if (omega_y == 0.0f) {
    element->R_matrix[2*BEAM_DOFS + 2] = 0;
    element->R_matrix[2*BEAM_DOFS + 3] = L;
    element->R_matrix[3*BEAM_DOFS + 2] = 0;
    element->R_matrix[3*BEAM_DOFS + 3] = 1;
  } else {
    if (omega_y_sqr > 0) {
      sinlike_func = sinh;
      coslike_func = cosh;
      sign = 1;
    } else {
      sinlike_func = sin;
      coslike_func = cos;
      sign = -1;
    }
    assert((sign == 1) | (sign == -1) && "Sign value has been corrupted");
    element->R_matrix[2*BEAM_DOFS + 2] = coslike_func(omega_y_L);
    element->R_matrix[2*BEAM_DOFS + 3] = (1/omega_y) * sinlike_func(omega_y_L);
    element->R_matrix[3*BEAM_DOFS + 2] = sign * omega_y * sinlike_func(omega_y_L);
    element->R_matrix[3*BEAM_DOFS + 3] = coslike_func(omega_y_L);
  }
}

void rmatrix_print(FILE *file, double mat[BEAM_DOFS*BEAM_DOFS]) {
  char *fmt_str = "%+0.6e";
  for (size_t j=0; j<BEAM_DOFS; j++) {
    for (size_t i=0; i<BEAM_DOFS; i++) {
      double val = mat[j*BEAM_DOFS + i];
      fprintf(file, fmt_str, val);
      if (i!=BEAM_DOFS-1) fprintf(file, ", ");
    }
    fprintf(file, "\n");
  }

  // double R = determinant(mat, BEAM_DOFS);
  // printf("|R| - 1 = %0.6e     ", R - 1);
  // if ((fabs(R - 1) > EPSILON)) printf("\n****** Non-unity determinant! ******");
  // printf("\n");
}

Element create_element(char *name, char **cursor) {
  Element result = {0};
  size_t name_len = strlen(name);
  memcpy(result.name, name, name_len >= ELENAME_MAX_LEN ? ELENAME_MAX_LEN-1 : name_len);

  Arena mem_arena = make_arena();
  str2dbl_hashmap* kv_pairs = NULL;

  while (**cursor == ' ' | **cursor == ':') (*cursor)++;

  char type[50] = {0};
  size_t i = 0;
  while (isalpha(**cursor)) {
    type[i++] = (char)toupper(**cursor);
    (*cursor)++;
  }

  if (strcmp(type, "DRIFT") == 0)           result.type = ELETYPE_DRIFT;
  else if (strcmp(type, "KICKER") == 0)     result.type = ELETYPE_DRIFT;
  else if (strcmp(type, "MONITOR") == 0)    result.type = ELETYPE_DRIFT;
  else if (strcmp(type, "MARKER") == 0)     result.type = ELETYPE_DRIFT;
  else if (strcmp(type, "SBEND") == 0)      result.type = ELETYPE_SBEND;
  else if (strcmp(type, "QUADRUPOLE") == 0) result.type = ELETYPE_QUAD;
  else if (strcmp(type, "SEXTUPOLE") == 0)  result.type = ELETYPE_SEXTUPOLE;
  else if (strcmp(type, "MULTIPOLE") == 0)  result.type = ELETYPE_MULTIPOLE;
  else if (strcmp(type, "RFCAVITY") == 0)   result.type = ELETYPE_CAVITY;
  else if (strcmp(type, "LINE") == 0) {
    *cursor -= 4; // Safe since we just extracted "LINE" from this, so this puts the cursor back there.
    arena_free(&mem_arena);
    return result;
  } else {
    fprintf(stderr, "ERROR: Unable to parse elements of type %s\n", type);
    printf("%s\n", *cursor);
    exit(1);
  }

  while (**cursor==' ' | **cursor==',') (*cursor)++;

  bool finding_key = true;
  bool finding_val = false;
  char *key = arena_alloc(&mem_arena, 100*sizeof(char));
  size_t key_index = 0;
  double val = 0.0;
  for (;;) {
    if (**cursor == ';') {
      break;
    }
    if (finding_key & !finding_val) {
      if (isvalididchar(**cursor)) {
        key[key_index] = **cursor;
        key_index++;
        (*cursor)++;
      } else {
        if (**cursor==' ' | **cursor=='=') {
          while (!(isdigit(**cursor) | (**cursor=='+') | (**cursor=='-')) & (**cursor!=';')) {
            (*cursor)++;
          }
          if (**cursor == ';') {
            fprintf(stderr, "Element definition terminated early");
            exit(1);
          }
          finding_key = false;
          finding_val = true;
        }
      }
    } else if (finding_val & !finding_key) {
      val = strtod(*cursor, cursor);
      shput(kv_pairs, key, val);
      finding_key = true;
      finding_val = false;
      key = arena_alloc(&mem_arena, 100*sizeof(char));
      key_index = 0;
      val = 0.0;
      while (**cursor==',' | **cursor==' ') {
        (*cursor)++;
      }
    }
  }

  if (kv_pairs == NULL) {
    shfree(kv_pairs);
    make_r_matrix(&result);
    arena_free(&mem_arena);
    return result;
  }

  switch (result.type) {
    case ELETYPE_DRIFT:
      result.as.drift.length = shget(kv_pairs, "L");
      break;
    case ELETYPE_QUAD:
      result.as.quad.length = shget(kv_pairs, "L");
      result.as.quad.K1     = shget(kv_pairs, "K1");
      break;
    case ELETYPE_SBEND:
      result.as.sbend.length = shget(kv_pairs, "L");
      result.as.sbend.K1     = shget(kv_pairs, "K1");
      result.as.sbend.angle  = shget(kv_pairs, "ANGLE");
      result.as.sbend.E1     = shget(kv_pairs, "E1");
      result.as.sbend.E2     = shget(kv_pairs, "E2");
      break;
    case ELETYPE_MULTIPOLE:
      result.as.multipole.length = shget(kv_pairs, "L");
      result.as.multipole.K1L    = shget(kv_pairs, "K1L");
      result.as.multipole.K2L    = shget(kv_pairs, "K2L");
      result.as.multipole.K3L    = shget(kv_pairs, "K3L");
      break;
    case ELETYPE_SEXTUPOLE:
      result.as.sextupole.length = shget(kv_pairs, "L");
      result.as.sextupole.K2     = shget(kv_pairs, "K2");
      break;
    case ELETYPE_CAVITY:
      result.as.cavity.length   = shget(kv_pairs, "L");
      result.as.cavity.voltage  = shget(kv_pairs, "VOLT");
      result.as.cavity.harmonic = shget(kv_pairs, "harm");
      result.as.cavity.lag      = shget(kv_pairs, "lag");
      break;
    case ELETYPE_OCTUPOLE:
      break;
  }

  Element new_result = {0};
  if (result.type == ELETYPE_MULTIPOLE) {
    double K1L = result.as.multipole.K1L;
    double K2L = result.as.multipole.K2L;
    double K3L = result.as.multipole.K3L;
    memcpy(new_result.name, result.name, strlen(result.name));
    if (K2L == 0.0 & K3L == 0.0) {
      // This is actually a quad
      new_result.type = ELETYPE_QUAD;
      new_result.as.quad.length = result.as.multipole.length;
      new_result.as.quad.K1 = K1L;
      result = new_result;
    } else if (K1L == 0.0 & K3L == 0.0) {
      // This is actually a sextupole
      new_result.type = ELETYPE_SEXTUPOLE;
      new_result.as.sextupole.length = result.as.multipole.length;
      new_result.as.sextupole.K2 = K2L;
      result = new_result;
    } else if (K1L == 0.0 & K2L == 0.0) {
      // This is actually an octupole
      new_result.type = ELETYPE_OCTUPOLE;
      new_result.as.octupole.length = result.as.multipole.length;
      new_result.as.octupole.K3 = K3L;
      result = new_result;
    } else if (K1L == 0.0 & K2L == 0.0 & K3L == 0.0) {
      // This is actually a drift
      new_result.type = ELETYPE_DRIFT;
      new_result.as.drift.length = result.as.multipole.length;
      result = new_result;
    }
  }

  make_r_matrix(&result);

  shfree(kv_pairs);
  arena_free(&mem_arena);

  return result;
}

bool isvalididchar(char c) {
  return isalnum(c) | (c == '_');
}

char *populate_element_library(char *cursor) {
  while (*cursor != '\0') {
    if (isalpha(*cursor)) {
      char *element_name = cursor;
      while (isvalididchar(*cursor)) {
        cursor++;
      }
      *cursor = '\0';
      cursor++;

      shput(element_library, element_name, arrput(element_list, 
                                                  create_element(element_name, &cursor)));
    } else {
      cursor++;
    }
    if (strncmp(cursor, "LINE", 4) == 0) {
      break;
    }
  }
  return cursor;
}

void element_print(Element element) {
  printf("%s: ", element.name);
  switch (element.type) {
    case ELETYPE_DRIFT:
      printf("Drift: L = %f\n", element.as.drift.length);
      break;
    case ELETYPE_SBEND:
      printf("SBend: L = %f, Angle = %f, K1 = %f, E1 = %f, E2 = %f\n", 
             element.as.sbend.length,
             element.as.sbend.angle,
             element.as.sbend.K1,
             element.as.sbend.E1,
             element.as.sbend.E2);
      break;
    case ELETYPE_QUAD:
      printf("Quad: L = %f, K1 = %f\n", 
             element.as.quad.length,
             element.as.quad.K1);
      break;
    case ELETYPE_MULTIPOLE:
      printf("Multipole: L = %f, K1L = %f, K2L = %f, K3L = %f\n", 
             element.as.multipole.length,
             element.as.multipole.K1L,
             element.as.multipole.K2L,
             element.as.multipole.K3L);
      break;
    case ELETYPE_SEXTUPOLE:
      printf("Sextupole: L = %f, K2 = %f\n", 
             element.as.sextupole.length,
             element.as.sextupole.K2);
      break;
    case ELETYPE_OCTUPOLE:
      printf("Sextupole: L = %f, K2 = %f\n", 
             element.as.octupole.length,
             element.as.octupole.K3);
      break;
    case ELETYPE_CAVITY:
      printf("Cavity; L = %f, V = %f, harmonic = %f, lag = %f\n",
             element.as.cavity.length,
             element.as.cavity.voltage,
             element.as.cavity.harmonic,
             element.as.cavity.lag);
      break;
  }
}

double calculate_line_angle(Element *line) {
  double total_angle = 0;
  for (size_t i=0; i<arrlenu(line); i++) {
    switch (line[i].type) {
      case ELETYPE_SBEND:
        total_angle += line[i].as.sbend.angle;
        break;
      case ELETYPE_QUAD:
      case ELETYPE_DRIFT:
      case ELETYPE_MULTIPOLE:
      case ELETYPE_SEXTUPOLE:
      case ELETYPE_OCTUPOLE:
      case ELETYPE_CAVITY:
        break;
    }
  }
  return total_angle;
}

double calculate_line_length(Element *line) {
  double total_length = 0;
  for (size_t i=0; i<arrlenu(line); i++) {
    double ele_len = element_length(line[i]);
    total_length += ele_len;
  }
  return total_length;
}

void create_line(char *cursor, Element **line) {
  advance_to_char(&cursor, '(');
  while (*cursor != ')') {
    while (!isalpha(*cursor) & (*cursor != ')')) {
      cursor++;
    }
    if (*cursor == ')') {
      break;
    }
    char *temp_cursor = cursor;
    while (isvalididchar(*temp_cursor)) {
      temp_cursor++;
    }
    *temp_cursor = '\0';
    Element ele = shget(element_library, cursor);
    arrput(*line, ele);
    cursor = temp_cursor + 1;
  }
}

double e_loss_per_turn(double I2, double gamma0) {
  return ERADIUS_TIMES_RESTMASS * I2 * pow(gamma0, 4);
}

bool matrix_multiply(double *mat1, double *mat2, double *result, size_t r1, size_t c1, size_t r2, size_t c2) {
  if (c1 != r2) {
    fprintf(stderr, "Matrix dimensions do not allow for multiplication");
    return false;
  }

  size_t c3 = c2;

  for (size_t row=0; row<r1; row++) {
    for (size_t col=0; col<c2; col++) {
      result[row*c3 + col] = 0;

      for (size_t k=0; k<c1; k++) {
        result[row*c3 + col] += mat1[row*c3 + k] * mat2[k*c3 + col];
      }
    }
  }

  return true;
}

void apply_matrix_n_times(double* result, double *matrix, size_t N) {
  memcpy(result, SIXBYSIX_IDENTITY, BEAM_DOFS*BEAM_DOFS*sizeof(double));

  for (size_t i=0; i<N; i++) {
    double temp_result[BEAM_DOFS*BEAM_DOFS] = {0};
    matrix_multiply(matrix, result, temp_result, 6, 6, 6, 6);
    memcpy(result, temp_result, BEAM_DOFS*BEAM_DOFS*sizeof(double));
  }
}

void get_line_matrix(double *matrix, Element *line) {
  memcpy(matrix, SIXBYSIX_IDENTITY, BEAM_DOFS*BEAM_DOFS*sizeof(double));

  for (size_t i=0; i<arrlenu(line); i++) {
    double temp_result[BEAM_DOFS*BEAM_DOFS] = {0};
    matrix_multiply(line[i].R_matrix, matrix, temp_result, BEAM_DOFS, BEAM_DOFS, BEAM_DOFS, BEAM_DOFS);
    memcpy(matrix, temp_result, BEAM_DOFS*BEAM_DOFS*sizeof(double));
  }
}

// Courtesy ChatGPT
// Function to get the minor matrix by removing row `r` and column `c`
static void getMinor(double* matrix, double* minor, int n, int r, int c) {
    int minorRow = 0, minorCol = 0;
    for (int i = 0; i < n; i++) {
        if (i == r) continue; // Skip the row `r`
        minorCol = 0;
        for (int j = 0; j < n; j++) {
            if (j == c) continue; // Skip the column `c`
            minor[minorRow * (n - 1) + minorCol] = matrix[i * n + j];
            minorCol++;
        }
        minorRow++;
    }
}

// Courtesy ChatGPT
// Recursive function to calculate the determinant
double determinant(double* matrix, int n) {
    // Base case: if matrix is 2x2, return determinant directly
    if (n == 2) return matrix[0] * matrix[3] - matrix[1] * matrix[2];
    
    double det = 0.0;
    double minor[(n - 1) * (n - 1)]; // Minor matrix will be (n-1)x(n-1)

    // Cofactor expansion along the first row (i = 0)
    for (int j = 0; j < n; j++) {
        // Get the minor matrix for element at (0, j)
        getMinor(matrix, minor, n, 0, j);

        // Compute cofactor (-1)^(i+j) * element * determinant of minor
        double cofactor = (j % 2 == 0 ? 1 : -1) * matrix[j] * determinant(minor, n - 1);

        // Add cofactor to the determinant
        det += cofactor;
    }

    return det;
}

