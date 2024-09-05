#ifndef ELEMENTS_H
#define ELEMENTS_H

#include <stdbool.h>
#include <stddef.h>

#define C 299792458.0f
#define ELECTRON_MASS 510998.9499961642f
#define BEAM_DOFS 6
#define ERADIUS_TIMES_RESTMASS 0.959976365e-9

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
  ELETYPE_DRIFT,
  ELETYPE_QUAD,
  ELETYPE_SBEND,
  ELETYPE_MULTIPOLE,
  ELETYPE_SEXTUPOLE,
  ELETYPE_OCTUPOLE,
  ELETYPE_CAVITY,
} EleType;

typedef struct {
  char name[ELENAME_MAX_LEN];
  EleType type;
  double R_matrix[BEAM_DOFS*BEAM_DOFS];
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
  char* key;
  double value;
} str2dbl_hashmap;

double synch_rad_integral_1(Element *line, int periodicity);
double synch_rad_integral_2(Element *line, int periodicity);
double synch_rad_integral_3(Element *line, int periodicity);
double element_length(Element element);
double bending_radius_of_element(Element element);
Element create_element(char *name, char **cursor);
void element_print(Element element);
void rmatrix_print(double mat[BEAM_DOFS*BEAM_DOFS]);
bool isvalididchar(char c);
char *populate_element_library(char *cursor);
double calculate_line_length(Element *line);
double calculate_line_angle(Element *line);
void create_line(char *cursor, Element **line);
double e_loss_per_turn(double I2, double gamma0);
bool matrix_multiply(double *mat1, double *mat2, double *result, size_t r1, size_t c1, size_t r2, size_t c2);

#endif // !ELEMENTS_H

