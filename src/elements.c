#include "sdm_lib.h"
#define _GNU_SOURCE
#include <math.h>
#include <errno.h>
#include <assert.h>
#include <string.h>
#include <strings.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <cblas.h>

#define STB_DS_IMPLEMENTATION
#include "stb_ds.h"

#include "elements.h"
#include "lib.h"

#define EPSILON 1e-9

static void calc_sbend_matrix(Element *element);
static void calc_quad_matrix(Element *element);

typedef struct {
  char *key;
  double value;
} ParsingVar;

ParsingVar *variables = NULL;

bool element_is_nonlinear(Element ele) {
  return ((ele.type == ELETYPE_OCTUPOLE) ||
          (ele.type == ELETYPE_SEXTUPOLE) ||
          (ele.type == ELETYPE_MULTIPOLE));
}

void track_thru(double *beam, size_t n_particles, Element element) {
  double temp_result[BEAM_DOFS] = {0};
  if (!element_is_nonlinear(element)) {
    matrix_multiply(element.R_matrix, beam, temp_result, BEAM_DOFS, BEAM_DOFS, BEAM_DOFS, n_particles);
    memcpy(beam, temp_result, BEAM_DOFS * sizeof(double));
  } else if (element.type == ELETYPE_SEXTUPOLE) {
    Element half_ele = (Element){
      .type = ELETYPE_DRIFT, 
      .as.drift = {.length = element_length(element) / 2.0}
    };
    make_r_matrix(&half_ele);

    matrix_multiply(half_ele.R_matrix, beam, temp_result, BEAM_DOFS, BEAM_DOFS, BEAM_DOFS, n_particles);
    memcpy(beam, temp_result, BEAM_DOFS * sizeof(double));

    double K2L = element.as.sextupole.K2 * element.as.sextupole.length;
    for (size_t i=0; i<n_particles; i++) {
      double x = beam[0*n_particles + i];
      double y = beam[2*n_particles + i];
      beam[1*n_particles + i] += -0.5 * K2L * (x*x - y*y);
      beam[3*n_particles + i] += K2L * x * y;
    }

    matrix_multiply(half_ele.R_matrix, beam, temp_result, BEAM_DOFS, BEAM_DOFS, BEAM_DOFS, n_particles);
    memcpy(beam, temp_result, BEAM_DOFS * sizeof(double));
  } else if (element.type == ELETYPE_OCTUPOLE) {
    Element half_ele = (Element){
      .type = ELETYPE_DRIFT, 
      .as.drift = {.length = element_length(element) / 2.0}
    };
    make_r_matrix(&half_ele);

    matrix_multiply(half_ele.R_matrix, beam, temp_result, BEAM_DOFS, BEAM_DOFS, BEAM_DOFS, n_particles);
    memcpy(beam, temp_result, BEAM_DOFS * sizeof(double));

    double K3L = element.as.octupole.K3 * element.as.octupole.length;
    for (size_t i=0; i<n_particles; i++) {
      double x = beam[0*n_particles + i];
      double y = beam[2*n_particles + i];
      beam[1*n_particles + i] += -(1.0/6.0) * K3L * (x*x*x - 3*x*y*y);
      beam[3*n_particles + i] += -(1.0/6.0) * K3L * (y*y*y - 3*x*x*y);
    }

    matrix_multiply(half_ele.R_matrix, beam, temp_result, BEAM_DOFS, BEAM_DOFS, BEAM_DOFS, n_particles);
    memcpy(beam, temp_result, BEAM_DOFS * sizeof(double));
  } else {
    assert(0 && "Can't track this element yet");
  }
}

void track(double *beam, size_t n_particles, Element *line, size_t n_elements) {
  for (size_t i=0; i<n_elements; i++) {
    track_thru(beam, n_particles, line[i]);
  }
}

double synch_rad_integral_2(Element *line) {
  double I_2 = 0;
  for (size_t i=0; i<arrlenu(line); i++) {
    double rho = bending_radius_of_element(line[i]);
    I_2 += element_length(line[i]) / pow(rho, 2);
  }

  return I_2;
}

double synch_rad_integral_3(Element *line) {
  double I_3 = 0;
  for (size_t i=0; i<arrlenu(line); i++) {
    double rho_abs = fabs(bending_radius_of_element(line[i]));
    I_3 += element_length(line[i]) / pow(rho_abs, 3);
  }

  return I_3;
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
  fprintf(stderr, "Reached unreachable code. Something is badly wrong. Exiting.\n");
  exit(1);
}

const char* get_element_type(EleType t) {
  switch (t) {
    case ELETYPE_DRIFT: return "DRIFT";
    case ELETYPE_QUAD: return "QUAD";
    case ELETYPE_SBEND: return "SBEND";
    case ELETYPE_MULTIPOLE: return "MULTIPOLE";
    case ELETYPE_SEXTUPOLE: return "SEXTUPOLE";
    case ELETYPE_OCTUPOLE: return "OCTUPOLE";
    case ELETYPE_CAVITY: return "CAVITY";
  }
  return "";
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

void make_r_matrix(Element *element) {
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

  element->eta_prop_matrix[0*3 + 0] = element->R_matrix[0*BEAM_DOFS + 0];
  element->eta_prop_matrix[0*3 + 1] = element->R_matrix[0*BEAM_DOFS + 1];
  element->eta_prop_matrix[0*3 + 2] = element->R_matrix[0*BEAM_DOFS + 5];
  element->eta_prop_matrix[1*3 + 0] = element->R_matrix[1*BEAM_DOFS + 0];
  element->eta_prop_matrix[1*3 + 1] = element->R_matrix[1*BEAM_DOFS + 1];
  element->eta_prop_matrix[1*3 + 2] = element->R_matrix[1*BEAM_DOFS + 5];
  element->eta_prop_matrix[2*3 + 0] = 0.0f;
  element->eta_prop_matrix[2*3 + 1] = 0.0f;
  element->eta_prop_matrix[2*3 + 2] = 1.0f;
}

static void calc_quad_matrix(Element *element) {
  assert(element->type == ELETYPE_QUAD && "This func should only be called for quads.");
  double omega = sqrt(fabs(element->as.quad.K1));
  double L = element->as.quad.length;
  double omega_L = omega * L;

  if (element->as.quad.K1 == 0) {
    element->R_matrix[0*BEAM_DOFS + 1] = L;
    element->R_matrix[2*BEAM_DOFS + 3] = L;
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
  double h = angle / L;

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
    assert(((sign == 1) || (sign == -1)) && "Sign value has been corrupted");
    element->R_matrix[0*BEAM_DOFS + 0] = coslike_func(omega_x_L);
    element->R_matrix[0*BEAM_DOFS + 1] = (1/omega_x) * sinlike_func(omega_x_L);
    element->R_matrix[1*BEAM_DOFS + 0] = sign * omega_x * sinlike_func(omega_x_L);
    element->R_matrix[1*BEAM_DOFS + 1] = coslike_func(omega_x_L);

    element->R_matrix[0*BEAM_DOFS + 5] = (h/omega_x_sqr) * (1 - coslike_func(omega_x_L));
    element->R_matrix[1*BEAM_DOFS + 5] = (h/omega_x) * sinlike_func(omega_x_L);

    element->R_matrix[4*BEAM_DOFS + 0] = element->R_matrix[1*BEAM_DOFS + 5];
    element->R_matrix[4*BEAM_DOFS + 1] = element->R_matrix[0*BEAM_DOFS + 5];

    element->R_matrix[4*BEAM_DOFS + 4] = 1;
    element->R_matrix[4*BEAM_DOFS + 5] = -pow(h,2)*(omega_x_L - sinlike_func(omega_x_L)) / pow(omega_x, 3);
    element->R_matrix[5*BEAM_DOFS + 4] = 0;
    element->R_matrix[5*BEAM_DOFS + 5] = 1;
  }

  if (omega_y == 0.0f) {
    element->R_matrix[2*BEAM_DOFS + 2] = 1;
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
    assert(((sign == 1) || (sign == -1)) && "Sign value has been corrupted");
    element->R_matrix[2*BEAM_DOFS + 2] = coslike_func(omega_y_L);
    element->R_matrix[2*BEAM_DOFS + 3] = (1/omega_y) * sinlike_func(omega_y_L);
    element->R_matrix[3*BEAM_DOFS + 2] = sign * omega_y * sinlike_func(omega_y_L);
    element->R_matrix[3*BEAM_DOFS + 3] = coslike_func(omega_y_L);
  }
}

void rmatrix_print(FILE *file, double mat[BEAM_DOFS*BEAM_DOFS]) {
  char *fmt_str = "%+0.6e";
  for (size_t j=0; j<BEAM_DOFS; j++) {
    fprintf(file, " ");
    for (size_t i=0; i<BEAM_DOFS; i++) {
      double val = mat[j*BEAM_DOFS + i];
      fprintf(file, fmt_str, val);
      if (i!=BEAM_DOFS-1) fprintf(file, "  ");
    }
    fprintf(file, "\n");
  }
}

Element create_element(sdm_arena_t *mem_arena, char *name, char **cursor) {
  Element result = {0};
  size_t name_len = strlen(name);
  memcpy(result.name, name, name_len >= ELENAME_MAX_LEN ? ELENAME_MAX_LEN-1 : name_len);

  str2dbl_hashmap* kv_pairs = NULL;

  while ((**cursor == ' ') | (**cursor == ':')) (*cursor)++;

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
    return result;
  } else {
    fprintf(stderr, "ERROR: Unable to parse elements of type %s\n", type);
    printf("%s\n", *cursor);
    exit(1);
  }

  while ((**cursor==' ') | (**cursor==',')) (*cursor)++;

  bool finding_key = true;
  bool finding_val = false;
  char *key = sdm_arena_alloc(mem_arena, 100*sizeof(char));
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
        if ((**cursor==' ') | (**cursor=='=')) {
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
      key = sdm_arena_alloc(mem_arena, 100*sizeof(char));
      key_index = 0;
      val = 0.0;
      while ((**cursor==',') | (**cursor==' ')) {
        (*cursor)++;
      }
    }
  }

  if (kv_pairs == NULL) {
    shfree(kv_pairs);
    make_r_matrix(&result);
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
    if ((K2L == 0.0) & (K3L == 0.0)) {
      // This is actually a quad
      new_result.type = ELETYPE_QUAD;
      new_result.as.quad.length = result.as.multipole.length;
      new_result.as.quad.K1 = K1L;
      result = new_result;
    } else if ((K1L == 0.0) & (K3L == 0.0)) {
      // This is actually a sextupole
      new_result.type = ELETYPE_SEXTUPOLE;
      new_result.as.sextupole.length = result.as.multipole.length;
      new_result.as.sextupole.K2 = K2L;
      result = new_result;
    } else if ((K1L == 0.0) & (K2L == 0.0)) {
      // This is actually an octupole
      new_result.type = ELETYPE_OCTUPOLE;
      new_result.as.octupole.length = result.as.multipole.length;
      new_result.as.octupole.K3 = K3L;
      result = new_result;
    } else if ((K1L == 0.0) & (K2L == 0.0) & (K3L == 0.0)) {
      // This is actually a drift
      new_result.type = ELETYPE_DRIFT;
      new_result.as.drift.length = result.as.multipole.length;
      result = new_result;
    }
  }

  make_r_matrix(&result);

  shfree(kv_pairs);

  return result;
}

char *populate_element_library(ElementLibrary **element_library, Element **element_list, char *cursor) {
  sdm_arena_t mem_arena = {0};
  sdm_arena_init(&mem_arena, SDM_ARENA_DEFAULT_CAP);

  while (*cursor != '\0') {
    if (*cursor == '!') {
      while (*cursor != '\n') cursor++;
      cursor++;
      continue;
    }
    if (isalpha(*cursor)) {
      char *element_name = cursor;
      while (isvalididchar(*cursor)) {
        cursor++;
      }
      *cursor = '\0';
      cursor++;

      if (strcmp(element_name, "RF_ON") == 0) {
        while (*cursor != '\n') cursor++;
        continue;
      }

      shput(*element_library, element_name, arrput(*element_list, 
                                                  create_element(&mem_arena, element_name, &cursor)));
    } else {
      cursor++;
    }
    if (strncmp(cursor, "LINE", 4) == 0) {
      break;
    }
  }

  sdm_arena_free(&mem_arena);
  return cursor;
}

void element_print(FILE *sink, Element element) {
  fprintf(sink, "%s: ", element.name);
  switch (element.type) {
    case ELETYPE_DRIFT:
      fprintf(sink, "Drift: L = %f\n", element.as.drift.length);
      break;
    case ELETYPE_SBEND:
      fprintf(sink, "SBend: L = %f, Angle = %f, K1 = %f, E1 = %f, E2 = %f\n", 
             element.as.sbend.length,
             element.as.sbend.angle,
             element.as.sbend.K1,
             element.as.sbend.E1,
             element.as.sbend.E2);
      break;
    case ELETYPE_QUAD:
      fprintf(sink, "Quad: L = %f, K1 = %f\n", 
             element.as.quad.length,
             element.as.quad.K1);
      break;
    case ELETYPE_MULTIPOLE:
      fprintf(sink, "Multipole: L = %f, K1L = %f, K2L = %f, K3L = %f\n", 
             element.as.multipole.length,
             element.as.multipole.K1L,
             element.as.multipole.K2L,
             element.as.multipole.K3L);
      break;
    case ELETYPE_SEXTUPOLE:
      fprintf(sink, "Sextupole: L = %f, K2 = %f\n", 
             element.as.sextupole.length,
             element.as.sextupole.K2);
      break;
    case ELETYPE_OCTUPOLE:
      fprintf(sink, "Sextupole: L = %f, K2 = %f\n", 
             element.as.octupole.length,
             element.as.octupole.K3);
      break;
    case ELETYPE_CAVITY:
      fprintf(sink, "Cavity; L = %f, V = %f, harmonic = %f, lag = %f\n",
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

void create_line(char *cursor, Element **line, ElementLibrary *element_library) {
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

double natural_emittance_x(double I2, double I4, double I5, double gamma0) {
  return C_Q * pow(gamma0, 2) * I5 / (I2 - I4);
}

double energy_spread(double I2, double I3, double I4, double gamma0) {
  return C_Q * pow(gamma0, 2) * I3 / (2*I2 + I4);
}

double cube(double val) {
  return pow(val, 3);
}

double sqr(double val) {
  return val*val;
}

double get_curlyH(Element element, double eta0, double etap0, double beta0, double alpha0) {
  assert(beta0 > 0.0 && "Beta of 0.0 or less does not make sense. Something has gone wrong.");
  double gamma0 = (1/beta0) * (1 + alpha0*alpha0);

  if (element.type != ELETYPE_SBEND || ((element.type == ELETYPE_SBEND) && element.as.sbend.angle == 0.0)) {
    return gamma0*eta0*eta0 + 2*alpha0*eta0*etap0 + beta0*etap0*etap0;
  }

  double K1 = element.as.sbend.K1;
  double L = element.as.sbend.length;
  double angle = element.as.sbend.angle;
  double h = fabs(angle / L);

  double k_sqr = pow(h, 2) + K1;
  double k = sqrt(fabs(k_sqr));

  double I5 = 0;
  double K = fabs(k_sqr);
  double psi = k * L;
  if (k_sqr > 0) {
    I5 =
        L*fabs(cube(h)) * (gamma0*sqr(eta0)+2e0*alpha0*eta0*etap0+beta0*sqr(etap0))
        - 2e0*pow(h, 4)/(pow(K, 3e0/2e0)) * (sqrt(K)*(alpha0*eta0+beta0*etap0)*(cos(psi)-1e0) + (gamma0*eta0+alpha0*etap0)*(psi-sin(psi)))
        + fabs(pow(h, 5))/(4e0*pow(K, 5e0/2e0)) * (2e0*alpha0*sqrt(K)*(4e0*cos(psi)-cos(2e0*psi)-3e0) + beta0*K*(2e0*psi-sin(2e0*psi)) +gamma0*(6e0*psi-8e0*sin(psi)+sin(2e0*psi)));
  } else {
    I5 = 
        L*fabs(cube(h)) * (gamma0*sqr(eta0)+2e0*alpha0*eta0*etap0+beta0*sqr(etap0))
        + 2e0*pow(h, 4)/(pow(K, 3e0/2e0)) * (sqrt(K)*(alpha0*eta0+beta0*etap0)*(cosh(psi)-1e0) + (gamma0*eta0+alpha0*etap0)*(psi-sinh(psi)))
        + fabs(pow(h, 5))/(4e0*pow(K, 5e0/2e0)) * (2e0*alpha0*sqrt(K)*(4e0*cosh(psi)-cosh(2e0*psi)-3e0) -beta0*K*(2e0*psi-sinh(2e0*psi)) +gamma0*(6e0*psi-8e0*sinh(psi)+sinh(2e0*psi)));
  }

  return I5 / (L*fabs(cube(h)));
}

void propagate_linear_optics(Element *line, double *line_matrix, LinOptsParams *lin_opt_params, double *I) {
  double S = 0.0;
  arrput(lin_opt_params->Ss, 0.0);

  get_line_matrix(line_matrix, line);

  const double R11 = line_matrix[0*BEAM_DOFS + 0];
  const double R12 = line_matrix[0*BEAM_DOFS + 1];
  const double R22 = line_matrix[1*BEAM_DOFS + 1];
  const double R33 = line_matrix[2*BEAM_DOFS + 2];
  const double R34 = line_matrix[2*BEAM_DOFS + 3];
  const double R44 = line_matrix[3*BEAM_DOFS + 3];
  const double R16 = line_matrix[0*BEAM_DOFS + 5];
  const double R56 = line_matrix[4*BEAM_DOFS + 5];

  const double eta_x  = R16 / (1 - R11);
  const double phi_x = acos((R11 + R22) / 2.0f);
  const double phi_y = acos((R33 + R44) / 2.0f);
  const double beta_x = fabs(R12 / sin(phi_x));
  const double beta_y = fabs(R34 / sin(phi_y));

  double eta_vec[3] = {eta_x, 0.0f, 1.0f};

  double twiss_mat[BEAM_DOFS*BEAM_DOFS] = {
    beta_x, 0,        0,      0,        0,     0,
    0,      1/beta_x, 0,      0,        0,     0,
    0,      0,        beta_y, 0,        0,     0,
    0,      0,        0,      1/beta_y, 0,     0,
    0,      0,        0,      0,        0,     0,
    0,      0,        0,      0,        0,     0,
  };

  arrput(lin_opt_params->element_beta_xs, beta_x);
  arrput(lin_opt_params->element_beta_ys, beta_y);
  arrput(lin_opt_params->element_etas, eta_vec[0]);
  arrput(lin_opt_params->element_etaps, eta_vec[1]);

  for (size_t i=0; i<arrlenu(line); i++) {
    arrput(lin_opt_params->element_curlyH, get_curlyH(line[i], eta_vec[0], eta_vec[1], twiss_mat[0*BEAM_DOFS + 0], -twiss_mat[0*BEAM_DOFS + 1]));

    S += element_length(line[i]);
    arrput(lin_opt_params->Ss, S);

    advance_twiss_matrix(twiss_mat, line[i]);

    arrput(lin_opt_params->element_beta_xs, twiss_mat[0*BEAM_DOFS + 0]);
    arrput(lin_opt_params->element_beta_ys, twiss_mat[2*BEAM_DOFS + 2]);

    double temp_eta_vec[3] = {0};
    matrix_multiply(line[i].eta_prop_matrix, eta_vec, temp_eta_vec, 3, 3, 3, 1);
    memcpy(eta_vec, temp_eta_vec, 3*sizeof(double));
    arrput(lin_opt_params->element_etas, eta_vec[0]);
    arrput(lin_opt_params->element_etaps, eta_vec[1]);
  }

  I[0] = R56;
  I[1] = synch_rad_integral_2(line);
  I[2] = synch_rad_integral_3(line);
  I[3] = 0.0;
  I[4] = 0.0;
  for (size_t i=0; i<arrlenu(line); i++) {
    if (line[i].type == ELETYPE_SBEND) {
      double angle = line[i].as.sbend.angle;
      double L = line[i].as.sbend.length;
      double K1 = line[i].as.sbend.K1;
      double h = angle / L;
      if (h == 0.0) continue;

      double (*sinlike_func)(double);
      double (*coslike_func)(double);
      double omega_sqr = pow(h, 2) + K1;
      double omega = sqrt(fabs(omega_sqr));
      double sign = 0.0;
      if (omega_sqr > 0.0) {
        sinlike_func = sin;
        coslike_func = cos;
        sign = 1.0;
      } else {
        sinlike_func = sinh;
        coslike_func = cosh;
        sign = -1.0;
      }
      double mean_eta = lin_opt_params->element_etas[i] * sinlike_func(omega*L) / (omega*L)
                + sign * lin_opt_params->element_etaps[i] * (1 - coslike_func(omega*L)) / (omega*omega*L)
                + sign* h * (omega*L - sinlike_func(omega*L))/(pow(omega,3)*L);

      I[3] += (mean_eta * h * L * (2*K1 + h*h));
      I[4] += (L * pow(fabs(h), 3) * lin_opt_params->element_curlyH[i]);
    }
  }
}

void advance_twiss_matrix(double *twiss_mat, Element element) {
  double intermed_twiss_mat[BEAM_DOFS*BEAM_DOFS] = {0};
  double input_twiss[BEAM_DOFS*BEAM_DOFS] = {0};
  memcpy(input_twiss, twiss_mat, BEAM_DOFS*BEAM_DOFS * sizeof(double));

  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,   BEAM_DOFS, BEAM_DOFS, BEAM_DOFS, 1.0, input_twiss,      BEAM_DOFS, element.R_matrix,   BEAM_DOFS, 0.0, intermed_twiss_mat, BEAM_DOFS);
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, BEAM_DOFS, BEAM_DOFS, BEAM_DOFS, 1.0, element.R_matrix, BEAM_DOFS, intermed_twiss_mat, BEAM_DOFS, 0.0, twiss_mat,          BEAM_DOFS);
}

bool matrix_multiply(double *mat1, double *mat2, double *result, size_t r1, size_t c1, size_t r2, size_t c2) {
  if (c1 != r2) {
    fprintf(stderr, "Matrix dimensions do not allow for multiplication");
    return false;
  }

  // C = alpha.A*B + beta.C
  //  void cblas_dgemm(CBLAS_LAYOUT layout, CBLAS_TRANSPOSE TransA,
  //           CBLAS_TRANSPOSE TransB, const CBLAS_INT M, const CBLAS_INT N,
  //           const CBLAS_INT K, const double alpha, const double *A,
  //           const CBLAS_INT lda, const double *B, const CBLAS_INT ldb,
  //           const double beta, double *C, const CBLAS_INT ldc);
  // enum CBLAS_LAYOUT layout = CblasRowMajor;
  // enum CBLAS_TRANSPOSE TransA = CblasNoTrans, TransB = CblasNoTrans;
  const int M = r1, N = c2, K = c1;
  const int lda = c1, ldb = c2;
  const double alpha = 1.0, beta = 0.0;

  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, alpha, mat1, lda, mat2, ldb, beta, result, ldb);

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

void generate_lattice_from_mad8_file(const char *filename, Element **line) {
  char *buffer = read_entire_file(filename);
  char *cursor = buffer;

  cursor = join_lines(cursor);

  Element *element_list = NULL;
  ElementLibrary *element_library = NULL;
  cursor = populate_element_library(&element_library, &element_list, cursor);

  create_line(cursor, line, element_library);

  arrfree(element_list);
  shfree(element_library);

  free(buffer);
}

typedef enum {
  LEXEME_TYPE_SYMBOL,
  LEXEME_TYPE_NUMBER,
  LEXEME_TYPE_ASSIGNMENT,
  LEXEME_TYPE_ADD,
  LEXEME_TYPE_MULT,
  LEXEME_TYPE_SUB,
  LEXEME_TYPE_DIV,
  LEXEME_TYPE_OPAREN,
  LEXEME_TYPE_CPAREN,
  LEXEME_TYPE_SEMICOLON,
  LEXEME_TYPE_COLON,
  LEXEME_TYPE_COMMA,
  LEXEME_TYPE_COUNT,
} LexemeType;

char *lexeme_strings[LEXEME_TYPE_COUNT] = {
  [LEXEME_TYPE_SYMBOL] =     "LEXEME_TYPE_SYMBOL",
  [LEXEME_TYPE_NUMBER] =     "LEXEME_TYPE_NUMBER",
  [LEXEME_TYPE_ASSIGNMENT] = "LEXEME_TYPE_ASSIGNMENT",
  [LEXEME_TYPE_ADD] =        "LEXEME_TYPE_ADD",
  [LEXEME_TYPE_MULT] =       "LEXEME_TYPE_MULT",
  [LEXEME_TYPE_SUB] =        "LEXEME_TYPE_SUB",
  [LEXEME_TYPE_DIV] =        "LEXEME_TYPE_DIV",
  [LEXEME_TYPE_OPAREN] =     "LEXEME_TYPE_OPAREN",
  [LEXEME_TYPE_CPAREN] =     "LEXEME_TYPE_CPAREN",
  [LEXEME_TYPE_SEMICOLON] =  "LEXEME_TYPE_SEMICOLON",
  [LEXEME_TYPE_COLON] =      "LEXEME_TYPE_COLON",
  [LEXEME_TYPE_COMMA] =      "LEXEME_TYPE_COMMA",
};

typedef struct {
  LexemeType type;
  sdm_string_view content;
} Lexeme;

size_t starts_with_float(const char *input) {
  char *endptr;

  errno = 0;
  strtod(input, &endptr);
  if (errno != 0) return 0;
  size_t retval = endptr - input;
  return retval;
}

typedef struct {
  size_t capacity;
  size_t length;
  Lexeme *data;
} LexemeArray;

typedef struct {
  size_t capacity;
  size_t length;
  double *data;
} DoubleArray;

int prec(Lexeme lex) {
  if (lex.type == LEXEME_TYPE_SUB) return 50;
  if (lex.type == LEXEME_TYPE_ADD) return 50;
  if (lex.type == LEXEME_TYPE_DIV) return 100;
  if (lex.type == LEXEME_TYPE_MULT) return 100;
  fprintf(stderr, "'prec' called with inappropriate lexeme");
  exit(1);
}

double evaluate_lexeme_string(Lexeme *lexemes, size_t *index) {
  DoubleArray return_stack = {0};
  SDM_ENSURE_ARRAY_MIN_CAP(return_stack, 1024);
  LexemeArray operator_stack = {0};
  SDM_ENSURE_ARRAY_MIN_CAP(operator_stack, 1024);

  Lexeme lexeme = lexemes[*index];
  
  bool done_parsing = false;
  while (((*index) < arrlenu(lexemes)) && !done_parsing) {
    if (lexeme.type == LEXEME_TYPE_SEMICOLON) {
      break;
    }
    switch (lexeme.type) {
      case LEXEME_TYPE_SEMICOLON:
      case LEXEME_TYPE_COMMA: {
        done_parsing = true;
      } break;
      case LEXEME_TYPE_COLON:
      case LEXEME_TYPE_COUNT:
      case LEXEME_TYPE_ASSIGNMENT: {
        fprintf(stderr, "Faulty expression\n");
        exit(1);
      } break;
      case LEXEME_TYPE_SYMBOL: {
        char *var_name = sdm_sv_to_cstr(lexeme.content);
        double val = shget(variables, var_name);
        SDM_ARRAY_PUSH(return_stack, val);
        free(var_name);
      } break;
      case LEXEME_TYPE_NUMBER: {
        SDM_ARRAY_PUSH(return_stack, strtod(lexeme.content.data, NULL));
      } break;
      case LEXEME_TYPE_MULT:
      case LEXEME_TYPE_DIV:
      case LEXEME_TYPE_ADD:
      case LEXEME_TYPE_SUB: {
        while ((operator_stack.length > 0) && (prec(operator_stack.data[operator_stack.length-1]) >= prec(lexeme))) {
          Lexeme popped = SDM_ARRAY_POP(operator_stack);
          assert(return_stack.length >= 2);
          double b = SDM_ARRAY_POP(return_stack);
          double a = SDM_ARRAY_POP(return_stack);
          if (popped.type == LEXEME_TYPE_ADD) SDM_ARRAY_PUSH(return_stack, a+b)
          else if (popped.type == LEXEME_TYPE_SUB) SDM_ARRAY_PUSH(return_stack, a-b)
          else if (popped.type == LEXEME_TYPE_MULT) SDM_ARRAY_PUSH(return_stack, a*b)
          else if (popped.type == LEXEME_TYPE_DIV) SDM_ARRAY_PUSH(return_stack, a/b)
          else {
            fprintf(stderr, "operator_stack corrupted\n");
            exit(1);
          }
        }
        SDM_ARRAY_PUSH(operator_stack, lexeme);
      } break;
      case LEXEME_TYPE_OPAREN:
      case LEXEME_TYPE_CPAREN: {
        fprintf(stderr, "Cannot handle parentheses yet\n");
        exit(1);
      } break;
    }
    (*index)++;
    lexeme = lexemes[*index];
  }

  while (operator_stack.length > 0) {
    Lexeme popped = SDM_ARRAY_POP(operator_stack);
    assert(return_stack.length >= 2);
    double b = SDM_ARRAY_POP(return_stack);
    double a = SDM_ARRAY_POP(return_stack);
    if (popped.type == LEXEME_TYPE_ADD) SDM_ARRAY_PUSH(return_stack, a+b)
    else if (popped.type == LEXEME_TYPE_SUB) SDM_ARRAY_PUSH(return_stack, a-b)
    else if (popped.type == LEXEME_TYPE_MULT) SDM_ARRAY_PUSH(return_stack, a*b)
    else if (popped.type == LEXEME_TYPE_DIV) SDM_ARRAY_PUSH(return_stack, a/b)
    else {
      fprintf(stderr, "operator_stack corrupted\n");
      exit(1);
    }
  }

  assert(return_stack.length == 1 && "Unused variables on return stack");
  return return_stack.data[0];
}

typedef struct {
  char *key;
  Element value;
} EleLibItem;

typedef struct {
  size_t capacity;
  size_t length;
  Element *data;
} Line;

typedef struct {
  char *key;
  Line value;
} LineItem;

void generate_lattice_from_tracy_file(const char *filename, Element **line) {
  EleLibItem *element_library = NULL;

  char *buffer = read_entire_file(filename);
  sdm_string_view file_contents = sdm_cstr_as_sv(buffer);

  Lexeme *lexemes = {0};
  sdm_sv_trim(&file_contents);
  while (file_contents.length > 0) {
    Lexeme lexeme = {0};

    size_t len = starts_with_float(file_contents.data);
    if (isalpha(file_contents.data[0])) {
      lexeme.type = LEXEME_TYPE_SYMBOL;
      size_t i = 0;
      while (isvalididchar(file_contents.data[i])) i++;
      lexeme.content = sdm_sized_str_as_sv(file_contents.data, i);
      file_contents.data += i;
      file_contents.length -= i;
    } else if (len > 0) {
      lexeme.type = LEXEME_TYPE_NUMBER;
      lexeme.content = sdm_sized_str_as_sv(file_contents.data, len);
      file_contents.data += len;
      file_contents.length -= len;
    } else if (file_contents.data[0] == '=') {
      lexeme.type = LEXEME_TYPE_ASSIGNMENT;
      lexeme.content = sdm_sized_str_as_sv(file_contents.data, 1);
      file_contents.data += 1;
      file_contents.length -= 1;
    } else if (file_contents.data[0] == ';') {
      lexeme.type = LEXEME_TYPE_SEMICOLON;
      lexeme.content = sdm_sized_str_as_sv(file_contents.data, 1);
      file_contents.data += 1;
      file_contents.length -= 1;
    } else if (file_contents.data[0] == '/') {
      lexeme.type = LEXEME_TYPE_DIV;
      lexeme.content = sdm_sized_str_as_sv(file_contents.data, 1);
      file_contents.data += 1;
      file_contents.length -= 1;
    } else if (file_contents.data[0] == ':') {
      lexeme.type = LEXEME_TYPE_COLON;
      lexeme.content = sdm_sized_str_as_sv(file_contents.data, 1);
      file_contents.data += 1;
      file_contents.length -= 1;
    } else if (file_contents.data[0] == ',') {
      lexeme.type = LEXEME_TYPE_COMMA;
      lexeme.content = sdm_sized_str_as_sv(file_contents.data, 1);
      file_contents.data += 1;
      file_contents.length -= 1;
    } else if (file_contents.data[0] == '*') {
      lexeme.type = LEXEME_TYPE_MULT;
      lexeme.content = sdm_sized_str_as_sv(file_contents.data, 1);
      file_contents.data += 1;
      file_contents.length -= 1;
    } else if (file_contents.data[0] == '-') {
      lexeme.type = LEXEME_TYPE_SUB;
      lexeme.content = sdm_sized_str_as_sv(file_contents.data, 1);
      file_contents.data += 1;
      file_contents.length -= 1;
    } else if (file_contents.data[0] == '(') {
      lexeme.type = LEXEME_TYPE_OPAREN;
      lexeme.content = sdm_sized_str_as_sv(file_contents.data, 1);
      file_contents.data += 1;
      file_contents.length -= 1;
    } else if (file_contents.data[0] == ')') {
      lexeme.type = LEXEME_TYPE_CPAREN;
      lexeme.content = sdm_sized_str_as_sv(file_contents.data, 1);
      file_contents.data += 1;
      file_contents.length -= 1;
    } else {
      fprintf(stderr, "Trying to parse an unknown character, %c\n", file_contents.data[0]);
      exit(1);
    }

    arrput(lexemes, lexeme);

    sdm_sv_trim(&file_contents);
  }

  LineItem *lines = NULL;

  size_t i = 0;
  while (i < arrlenu(lexemes)) {
    Lexeme lexeme = lexemes[i];

    if (lexeme.type != LEXEME_TYPE_SYMBOL) {
      fprintf(stderr, "Expected lexeme type %s but got %s \""SDM_SV_F"\"\n", lexeme_strings[LEXEME_TYPE_SYMBOL], lexeme_strings[lexeme.type], SDM_SV_Vals(lexeme.content));
      exit(1);
    }
    char *name = sdm_sv_to_cstr(sdm_sized_str_as_sv(lexeme.content.data, lexeme.content.length));

    if (strncasecmp(name, "USE", 3) == 0) {
      i++;
      lexeme = lexemes[i];
      assert(lexeme.type == LEXEME_TYPE_COLON);
      i++;
      lexeme = lexemes[i];
      assert(lexeme.type == LEXEME_TYPE_SYMBOL);
      char *line_name_to_use = sdm_sv_to_cstr(lexeme.content);
      Line line_to_use = shget(lines, line_name_to_use);
      for (size_t i=0; i<line_to_use.length; i++) {
        arrput((*line), line_to_use.data[i]);
      }
      i++;
      lexeme = lexemes[i];
      assert(lexeme.type == LEXEME_TYPE_SEMICOLON);
      i++;
      lexeme = lexemes[i];
      continue;
    }

    i++;
    lexeme = lexemes[i];

    if (lexeme.type == LEXEME_TYPE_ASSIGNMENT) {
      // We are assigning a value to a variable
      i++;
      lexeme = lexemes[i];
      double val = evaluate_lexeme_string(lexemes, &i);
      shput(variables, name, val);
    } else if (lexeme.type == LEXEME_TYPE_COLON) {
      // We are defining an element
      ParsingVar *ele_defns = NULL;
      i++;
      lexeme = lexemes[i];
      EleType ele_type = {0};
      if (strncasecmp(lexeme.content.data, "cavity", 6) == 0) ele_type = ELETYPE_CAVITY;
      else if (strncasecmp(lexeme.content.data, "drift", 5) == 0) ele_type = ELETYPE_DRIFT;
      else if (strncasecmp(lexeme.content.data, "quadrupole", 10) == 0) ele_type = ELETYPE_QUAD;
      else if (strncasecmp(lexeme.content.data, "bending", 7) == 0) ele_type = ELETYPE_SBEND;
      else if (strncasecmp(lexeme.content.data, "sextupole", 9) == 0) ele_type = ELETYPE_SEXTUPOLE;
      else if (strncasecmp(lexeme.content.data, "octupole", 8) == 0) ele_type = ELETYPE_OCTUPOLE;
      else if (strncasecmp(lexeme.content.data, "marker", 6) == 0) ele_type = ELETYPE_DRIFT;
      else if (strncasecmp(lexeme.content.data, "line", 4) == 0) {
        // Dealing with a line
        bool parse_backwards = false;
        Line newline = {0};
        SDM_ENSURE_ARRAY_MIN_CAP(newline, 256);
        i++;
        lexeme = lexemes[i];
        assert(lexeme.type == LEXEME_TYPE_ASSIGNMENT);
        i++;
        lexeme = lexemes[i];
        assert(lexeme.type = LEXEME_TYPE_OPAREN);
        i++;
        lexeme = lexemes[i];
        while (i < arrlenu(lexemes)) {
          if (lexeme.type == LEXEME_TYPE_COMMA) {
            parse_backwards = false;
            i++;
            lexeme= lexemes[i];
            continue;
          }
          else if (lexeme.type == LEXEME_TYPE_SYMBOL) {
            char *symbol_val = sdm_sv_to_cstr(lexeme.content);
            if (shgeti(element_library, symbol_val) >= 0) {
              SDM_ARRAY_PUSH(newline, shget(element_library, symbol_val));
            } else if (shgeti(lines, symbol_val) >= 0) {
              Line line2insert = shget(lines, symbol_val);
              if (parse_backwards) {
                for (int i=line2insert.length-1; i>=0; i--) {
                  SDM_ARRAY_PUSH(newline, line2insert.data[i]);
                }
              } else {
                for (size_t i=0; i<line2insert.length; i++) {
                  SDM_ARRAY_PUSH(newline, line2insert.data[i]);
                }
              }
            } else {
              fprintf(stderr, "%s is neither in element_library or lines\n", symbol_val);
              exit(1);
            }
            free(symbol_val);
          } else if (lexeme.type == LEXEME_TYPE_SUB) {
            parse_backwards = true;
          } else if (lexeme.type == LEXEME_TYPE_CPAREN) {
            i++;
            lexeme = lexemes[i];
            break;
          }
          i++;
          lexeme = lexemes[i];
        }
        shput(lines, name, newline);
        i++;
        lexeme = lexemes[i];
        continue;
      } else {
        fprintf(stderr, "Unknown element type: '"SDM_SV_F"'\n", SDM_SV_Vals(lexeme.content));
        exit(1);
      }
      i++;
      lexeme = lexemes[i];
      while (lexeme.type != LEXEME_TYPE_SEMICOLON) {
        while (lexeme.type != LEXEME_TYPE_SYMBOL) {
          i++;
          lexeme = lexemes[i];
        }
        char *varname = sdm_sv_to_cstr(sdm_sized_str_as_sv(lexeme.content.data, lexeme.content.length));
        i++;
        lexeme = lexemes[i];
        assert(lexeme.type == LEXEME_TYPE_ASSIGNMENT);
        i++;
        lexeme = lexemes[i];
        double val = evaluate_lexeme_string(lexemes, &i);
        shput(ele_defns, varname, val);
        lexeme = lexemes[i];
      }
      Element ele = {0};
      memcpy(ele.name, name, strlen(name)<64 ? strlen(name) : 63);
      switch (ele_type) {
        case ELETYPE_QUAD: {
          ele.type = ELETYPE_QUAD;
          ele.as.quad.length = shget(ele_defns, "L");
          ele.as.quad.K1 = shget(ele_defns, "B_2");
        } break;
        case ELETYPE_SBEND: {
          ele.type = ELETYPE_SBEND;
          ele.as.sbend.length = shget(ele_defns, "L");
          ele.as.sbend.K1 = shget(ele_defns, "B_2");
          ele.as.sbend.angle = degrees_to_radians(shget(ele_defns, "Phi"));
        } break;
        case ELETYPE_SEXTUPOLE: {
          ele.type = ELETYPE_SEXTUPOLE;
          ele.as.sextupole.length = shget(ele_defns, "L");
          ele.as.sextupole.K2 = shget(ele_defns, "B_3");
        } break;
        case ELETYPE_OCTUPOLE: {
          ele.type = ELETYPE_OCTUPOLE;
          ele.as.octupole.length = shget(ele_defns, "L");
          ele.as.octupole.K3 = shget(ele_defns, "B_4");
        } break;
        case ELETYPE_MULTIPOLE: {
          fprintf(stderr, "%s not implmented\n", get_element_type(ele_type));
          exit(1);
        } break;
        case ELETYPE_DRIFT: {
          ele.type = ELETYPE_DRIFT;
          ele.as.drift.length = shget(ele_defns, "L");
        } break;
        case ELETYPE_CAVITY: {
          ele.type = ELETYPE_CAVITY;
          ele.as.cavity.length   = shget(ele_defns, "L");
          ele.as.cavity.voltage = shget(ele_defns, "Voltage");
          ele.as.cavity.harmonic = shget(ele_defns, "HarNum");
        } break;
      }
      make_r_matrix(&ele);
      shput(element_library, name, ele);
    } else {
      fprintf(stderr, "Malformed input file\n");
      exit(1);
    }

    i++;
  }

  arrfree(lexemes);
  free(buffer);
}

