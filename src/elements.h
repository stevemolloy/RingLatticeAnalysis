#ifndef ELEMENTS_H
#define ELEMENTS_H

#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>

#include "sdm_lib.h"

#define C 299792458.0f
#define ELECTRON_MASS 510998.9499961642f
#define BEAM_DOFS 6
#define ERADIUS_TIMES_RESTMASS 0.959976365e-9
#define C_Q 3.83193864121903e-13

#define SIXBYSIX_IDENTITY (double[]){ \
    1, 0, 0, 0, 0, 0, \
    0, 1, 0, 0, 0, 0, \
    0, 0, 1, 0, 0, 0, \
    0, 0, 0, 1, 0, 0, \
    0, 0, 0, 0, 1, 0, \
    0, 0, 0, 0, 0, 1, \
  }

#define ELENAME_MAX_LEN 64

typedef struct {
  double length;
} Drift;

typedef struct {
  double length;
  double K1;
} Quad;

// SBEND,L=0.20424000,  ANGLE = 0.00388598, K1 = -1.33342000  E1= 0.00000000,  E2= 0.00000000 ; 
typedef struct {
  double length;
  double angle;
  double K1;
  double E1;
  double E2;
} Sbend;

// Multipole, K3L= -38933.64000000; 
typedef struct {
  double length;
  double K1L;
  double K2L;
  double K3L;
} Multipole;

//Sextupole,  L= 0.10000000, K2= -272.67400000 ;
typedef struct {
  double length;
  double K2;
} Sextupole;

typedef struct {
  double length;
  double K3;
} Octupole;

// CAV : RFCavity, L=0.37800000,VOLT=0.09000000, harm=176, lag=0.0; 
typedef struct {
  double length;
  double voltage;
  double harmonic;
  double lag;
} Cavity;

typedef enum {
  ELETYPE_DRIFT = 0,
  ELETYPE_QUAD,
  ELETYPE_SBEND,
  ELETYPE_MULTIPOLE,
  ELETYPE_SEXTUPOLE,
  ELETYPE_OCTUPOLE,
  ELETYPE_CAVITY,
} EleType;

const char* get_element_type(EleType t);

typedef struct {
  char name[ELENAME_MAX_LEN];
  EleType type;
  double R_matrix[BEAM_DOFS*BEAM_DOFS];
  double eta_prop_matrix[9];
  union {
    Drift drift;
    Quad quad;
    Sbend sbend;
    Multipole multipole;
    Sextupole sextupole;
    Octupole octupole;
    Cavity cavity;
  } as;
} Element;

typedef struct {
  char* key;
  Element value;
} ElementLibrary;

typedef struct {
  double *Ss;
  double *element_beta_xs;
  double *element_beta_ys;
  double *element_etas;
  double *element_etaps;
  double *element_curlyH;
} LinOptsParams;

typedef enum {
  TOKEN_TYPE_SYMBOL,
  TOKEN_TYPE_NUMBER,
  TOKEN_TYPE_ASSIGNMENT,
  TOKEN_TYPE_ADD,
  TOKEN_TYPE_MULT,
  TOKEN_TYPE_SUB,
  TOKEN_TYPE_DIV,
  TOKEN_TYPE_OPAREN,
  TOKEN_TYPE_CPAREN,
  TOKEN_TYPE_SEMICOLON,
  TOKEN_TYPE_COLON,
  TOKEN_TYPE_COMMA,
  TOKEN_TYPE_COUNT,
} TokenType;

typedef struct {
  TokenType type;
  sdm_string_view content;
} Token;

typedef struct {
  size_t capacity;
  size_t length;
  Token *data;
} TokenArray;

typedef struct {
  size_t capacity;
  size_t length;
  double *data;
} DoubleArray;

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

void track_thru(double *beam, size_t n_particles, Element element);
void track(double *beam, size_t n_particles, Element *line, size_t n_elements);

double synch_rad_integral_2(Element *line);
double synch_rad_integral_3(Element *line);
double synch_rad_integral_4(Element *line, int periodicity, double *element_etas);
double synch_rad_integral_5(Element *line, int periodicity, double *element_curlyH);
double e_loss_per_turn(double I2, double gamma0);
double natural_emittance_x(double I2, double I4, double I5, double gamma0);
double energy_spread(double I2, double I3, double I4, double gamma0);
double get_curlyH(Element element, double eta, double etap, double beta, double alpha);
void propagate_linear_optics(Element *line, double *total_matrix, LinOptsParams *lin_opt_params, double *I_synch);

void generate_lattice_from_mad8_file(const char *filename, Element **line);
void generate_lattice_from_tracy_file(const char *filename, Element **line);
bool tokenise_tracy_file(sdm_string_view *file_contents, Token **tokens);
void create_line(char *cursor, Element **line, ElementLibrary *element_library);
char *populate_element_library(ElementLibrary **element_library, Element **element_list, char *cursor);
void get_line_matrix(double *matrix, Element *line);
Element create_element(sdm_arena_t *mem_arena, char *name, char **cursor);
void element_print(FILE *sink, Element element);
double element_length(Element element);
double bending_radius_of_element(Element element);
double calculate_line_length(Element *line);
double calculate_line_angle(Element *line);
void make_r_matrix(Element *element);

void rmatrix_print(FILE *file, double mat[BEAM_DOFS*BEAM_DOFS]);
bool matrix_multiply(double *mat1, double *mat2, double *result, size_t r1, size_t c1, size_t r2, size_t c2);
void advance_twiss_matrix(double *twiss_mat, Element element);
void apply_matrix_n_times(double* result, double *matrix, size_t N);

#endif // !ELEMENTS_H

