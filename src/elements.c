#define _GNU_SOURCE
#include <math.h>
#include <assert.h>
#include <string.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>

#include "stb_ds.h"

#include "elements.h"
#include "lib.h"

extern Element* element_list;
extern ElementLibrary *element_library;

double synch_rad_integral_1(Element *line, size_t periodicity) {
  (void)line;
  (void)periodicity;
  return 0;
}

double synch_rad_integral_2(Element *line, size_t periodicity) {
  double I_2 = 0;
  for (size_t i=0; i<arrlenu(line); i++) {
    double rho = bending_radius_of_element(line[i]);
    I_2 += element_length(line[i]) / pow(rho, 2);
  }

  return I_2 * periodicity;
}

double synch_rad_integral_3(Element *line, size_t periodicity) {
  double I_3 = 0;
  for (size_t i=0; i<arrlenu(line); i++) {
    double rho_abs = fabs(bending_radius_of_element(line[i]));
    I_3 += element_length(line[i]) / pow(rho_abs, 3);
  }

  return I_3 * periodicity;
}

double element_length(Element element) {
  double length;
  switch (element.type) {
    case ELETYPE_SBEND:
      length = element.as.sbend.length;
      break;
    case ELETYPE_DRIFT:
      length = element.as.drift.length;
      break;
    case ELETYPE_QUAD:
      length = element.as.quad.length;
      break;
    case ELETYPE_SEXTUPOLE:
      length = element.as.sextupole.length;
      break;
    case ELETYPE_MULTIPOLE:
      length = element.as.multipole.length;
      break;
  }
  return length;
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
    case ELETYPE_MULTIPOLE:
      break;
  }
  return rho;
}

void make_r_matrix(Element *element) {
  double ele_length = element_length(*element);
  double omega = 0;

  memset(element->R_matrix, 0, sizeof(element->R_matrix));
  for (size_t i=0; i<BEAM_DOFS; i++) {
    element->R_matrix[i*BEAM_DOFS + i] = 1.0;
  }

  switch (element->type) {
    case (ELETYPE_DRIFT):
    case (ELETYPE_MULTIPOLE):
    case (ELETYPE_SEXTUPOLE):
      element->R_matrix[0*BEAM_DOFS + 1] = ele_length;
      element->R_matrix[2*BEAM_DOFS + 3] = ele_length;
      break;
    case (ELETYPE_QUAD):
      if (element->as.quad.K1 == 0) {
        element->R_matrix[0*BEAM_DOFS + 1] = ele_length;
        element->R_matrix[0*BEAM_DOFS + 1] = ele_length;
        break;
      }

      omega = sqrt(fabs(element->as.quad.K1));
      if (element->as.quad.K1 > 0) {
        // Focusing 2x2
        element->R_matrix[0*BEAM_DOFS + 0] = cos(omega * ele_length);
        element->R_matrix[0*BEAM_DOFS + 1] = sin(omega * ele_length) / omega;
        element->R_matrix[1*BEAM_DOFS + 0] = sin(omega * ele_length) * (-omega);
        element->R_matrix[1*BEAM_DOFS + 1] = cos(omega * ele_length);
        // Defocusing 2x2
        element->R_matrix[2*BEAM_DOFS + 2] = cosh(omega * ele_length);
        element->R_matrix[2*BEAM_DOFS + 3] = sinh(omega * ele_length) / omega;
        element->R_matrix[3*BEAM_DOFS + 2] = sinh(omega * ele_length) * omega;
        element->R_matrix[3*BEAM_DOFS + 3] = cosh(omega * ele_length);
      } else if (element->as.quad.K1 < 0) {
        // Defocusing 2x2
        element->R_matrix[2*BEAM_DOFS + 2] = cos(omega * ele_length);
        element->R_matrix[2*BEAM_DOFS + 3] = sin(omega * ele_length) / omega;
        element->R_matrix[3*BEAM_DOFS + 2] = sin(omega * ele_length) * (-omega);
        element->R_matrix[3*BEAM_DOFS + 3] = cos(omega * ele_length);
        // Focusing 2x2
        element->R_matrix[0*BEAM_DOFS + 0] = cosh(omega * ele_length);
        element->R_matrix[0*BEAM_DOFS + 1] = sinh(omega * ele_length) / omega;
        element->R_matrix[1*BEAM_DOFS + 0] = sinh(omega * ele_length) * omega;
        element->R_matrix[1*BEAM_DOFS + 1] = cosh(omega * ele_length);
      }

      break;
    case (ELETYPE_SBEND):
      calc_sbend_matrix(element);
      break;
  }
}

void calc_sbend_matrix(Element *element) {
  assert(element->type == ELETYPE_SBEND && "This func should only be called for sbends.");
  float K = element->as.sbend.K1;
  float L = element->as.sbend.length;
  float h = element->as.sbend.angle / L;

  float k_x_sqr = (1 - K) * h*h;
  float k_x = sqrt(fabs(k_x_sqr));
  if (k_x == 0.0f) {
    k_x = DBL_MIN;
    printf("k_x = %f\n", k_x);
  }
  float kxL = k_x * L;

  float k_y_sqr = -K * h*h;
  float k_y = sqrt(fabs(k_y_sqr));
  if (k_y == 0.0f) k_y = DBL_MIN;
  float kyL = k_y * L;

  double (*sinlike_func)(double);
  double (*coslike_func)(double);

  if (k_x_sqr > 0) {
    sinlike_func = &sin;
    coslike_func = &cos;
  } else {
    sinlike_func = &sinh;
    coslike_func = &cosh;
  }
  if (k_x == 0.0f) {
    element->R_matrix[0*BEAM_DOFS + 0] = 1;
    element->R_matrix[0*BEAM_DOFS + 1] = L;
    element->R_matrix[1*BEAM_DOFS + 0] = 0;
    element->R_matrix[1*BEAM_DOFS + 1] = 1;

    element->R_matrix[4*BEAM_DOFS + 4] = 1;
    element->R_matrix[4*BEAM_DOFS + 5] = 0;
    element->R_matrix[5*BEAM_DOFS + 4] = 0;
    element->R_matrix[5*BEAM_DOFS + 5] = 1;

    element->R_matrix[0*BEAM_DOFS + 5] = 0;
    element->R_matrix[1*BEAM_DOFS + 5] = h * L;
  } else {
    element->R_matrix[0*BEAM_DOFS + 0] =           coslike_func(kxL);
    element->R_matrix[0*BEAM_DOFS + 1] = (1/k_x) * sinlike_func(kxL);
    element->R_matrix[1*BEAM_DOFS + 0] =    -k_x * sinlike_func(kxL);
    element->R_matrix[1*BEAM_DOFS + 1] =           coslike_func(kxL);

    element->R_matrix[4*BEAM_DOFS + 4] = 1;
    element->R_matrix[4*BEAM_DOFS + 5] = (h*h/(k_x*k_x*k_x)) * (kxL - sinlike_func(kxL));
    element->R_matrix[5*BEAM_DOFS + 4] = 0;
    element->R_matrix[5*BEAM_DOFS + 5] = 1;

    element->R_matrix[0*BEAM_DOFS + 5] = (h/(k_x*k_x)) * (1 - coslike_func(kxL));
    element->R_matrix[1*BEAM_DOFS + 5] = (h/k_x) * sinlike_func(kxL);
  }

  if (k_y_sqr > 0) {
    sinlike_func = sin;
    coslike_func = cos;
  } else {
    sinlike_func = sinh;
    coslike_func = cosh;
  }
  if (k_y == 0.0f) {
    element->R_matrix[2*BEAM_DOFS + 2] = 0;
    element->R_matrix[2*BEAM_DOFS + 3] = L;
    element->R_matrix[3*BEAM_DOFS + 2] = 0;
    element->R_matrix[3*BEAM_DOFS + 3] = 1;

    element->R_matrix[4*BEAM_DOFS + 0] = L;
    element->R_matrix[4*BEAM_DOFS + 1] = 0;
  } else {
    element->R_matrix[2*BEAM_DOFS + 2] = coslike_func(kyL);
    element->R_matrix[2*BEAM_DOFS + 3] = (1/k_y) * sinlike_func(kyL);
    element->R_matrix[3*BEAM_DOFS + 2] = -k_y * sinlike_func(kyL);
    element->R_matrix[3*BEAM_DOFS + 3] = coslike_func(kyL);

    element->R_matrix[4*BEAM_DOFS + 0] = (1 / k_y) * sinlike_func(kyL);
    element->R_matrix[4*BEAM_DOFS + 1] = (h/(k_x*k_x)) * (1 - coslike_func(kxL));
  }
}

void rmatrix_print(double mat[BEAM_DOFS*BEAM_DOFS]) {
  for (size_t j=0; j<BEAM_DOFS; j++) {
    for (size_t i=0; i<BEAM_DOFS; i++) {
      double val = mat[j*BEAM_DOFS + i];
      (i==BEAM_DOFS-1) ? printf("%0.3f", val) : printf("%0.3f, ", val);
    }
    printf("\n");
  }
}

Element create_element(char **cursor) {
  char type[10] = {0};
  Element result = {0};

  while (**cursor == ' ' | **cursor == ':') (*cursor)++;

  size_t i = 0;
  while (isalpha(**cursor)) {
    type[i++] = toupper(**cursor);
    (*cursor)++;
  }

  if (strcmp(type, "DRIFT") == 0 | strcmp(type, "KICKER") == 0) {
    result.type = ELETYPE_DRIFT;
    char *temp_cursor = *cursor;
    while (*temp_cursor != ';') {
      temp_cursor++;
    }
    *temp_cursor = '\0';
    char *l_string = strcasestr(*cursor, "l");
    result.as.drift.length = assigned_double_from_string(l_string);
    while (**cursor != '\0') (*cursor)++;
    (*cursor)++;
  } else if (strcmp(type, "MARKER") == 0) {
    result.type = ELETYPE_DRIFT;
    result.as.drift.length = 0.0;
  } else if (strcmp(type, "SBEND") == 0) {
    result.type = ELETYPE_SBEND;
    char *temp_cursor = *cursor;
    while (*temp_cursor != ';') {
      temp_cursor++;
    }
    *temp_cursor = '\0';
    char *l_string = strcasestr(*cursor, "l");
    char *angle_string = strcasestr(*cursor, "angle");
    char *k1_string = strcasestr(*cursor, "k1");
    char *e1_string = strcasestr(*cursor, "e1");
    char *e2_string = strcasestr(*cursor, "e2");
    result.as.sbend.length = assigned_double_from_string(l_string);
    result.as.sbend.angle = assigned_double_from_string(angle_string);
    result.as.sbend.K1 = assigned_double_from_string(k1_string);
    result.as.sbend.E1 = assigned_double_from_string(e1_string);
    result.as.sbend.E2 = assigned_double_from_string(e2_string);
    while (**cursor != '\0') (*cursor)++;
    (*cursor)++;
  } else if (strcmp(type, "QUADRUPOLE") == 0) {
    result.type = ELETYPE_QUAD;
    char *temp_cursor = *cursor;
    while (*temp_cursor != ';') {
      temp_cursor++;
    }
    *temp_cursor = '\0';
    char *l_string = strcasestr(*cursor, "l");
    char *k1_string = strcasestr(*cursor, "k1");
    result.as.quad.length = assigned_double_from_string(l_string);
    result.as.quad.K1 = assigned_double_from_string(k1_string);
    while (**cursor != '\0') (*cursor)++;
    (*cursor)++;
  } else if (strcmp(type, "SEXTUPOLE") == 0) {
    result.type = ELETYPE_SEXTUPOLE;
    char *temp_cursor = *cursor;
    while (*temp_cursor != ';') {
      temp_cursor++;
    }
    *temp_cursor = '\0';
    char *l_string = strcasestr(*cursor, "l");
    char *k2_string = strcasestr(*cursor, "k2");
    result.as.sextupole.length = assigned_double_from_string(l_string);
    result.as.sextupole.K2 = assigned_double_from_string(k2_string);
    while (**cursor != '\0') (*cursor)++;
    (*cursor)++;
  } else if (strcmp(type, "MULTIPOLE") == 0) {
    result.type = ELETYPE_MULTIPOLE;
    char *temp_cursor = *cursor;
    while (*temp_cursor != ';') {
      temp_cursor++;
    }
    *temp_cursor = '\0';
    char *l_string = strcasestr(*cursor, "l");
    char *k3l_string = strcasestr(*cursor, "k3l");
    result.as.multipole.length = assigned_double_from_string(l_string);
    result.as.multipole.K3L = assigned_double_from_string(k3l_string);
    while (**cursor != '\0') (*cursor)++;
    (*cursor)++;
  } else if (strcmp(type, "LINE") == 0) {
    *cursor -= 4; // Safe since we just extracted "LINE" from this, so this puts the cursor back there.
    return result;
  } else {
    fprintf(stderr, "ERROR: Unable to parse elements of type %s\n", type);
    printf("%s\n", *cursor);
    exit(1);
  }

  make_r_matrix(&result);
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

      shput(element_library, element_name, arrput(element_list, create_element(&cursor)));
    } else {
      cursor++;
    }
    if (strncmp(cursor, "LINE", 4) == 0) {
      break;
    }
  }
  return cursor;
}

double assigned_double_from_string(char *string) {
  if (string) {
    while (*string != '=') string++;
    while ((*string != '-') & !isdigit(*string)) string++;
    return strtof(string, NULL);
  }
  return 0;
}

void element_print(Element element) {
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
      printf("Multipole: L = %f, K3L = %f\n", element.as.multipole.length, element.as.multipole.K3L);
      break;
    case ELETYPE_SEXTUPOLE:
      printf("Sextupole: L = %f, K2 = %f\n", 
             element.as.sextupole.length,
             element.as.sextupole.K2);
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

inline float e_loss_per_turn(float I2, float gamma0) {
  return ERADIUS_TIMES_RESTMASS * I2 * pow(gamma0, 4);
}

