#ifndef ELEMENTS_H
#define ELEMENTS_H

#include <stdbool.h>
#include <stddef.h>

#define C 299792458.0f
#define ELECTRON_MASS 510998.9499961642f
#define BEAM_DOFS 6
#define ERADIUS_TIMES_RESTMASS 0.959976365e-9

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
  double K3L;
} Multipole;

//Sextupole,  L= 0.10000000, K2= -272.67400000 ;
typedef struct {
  double length;
  double K2;
} Sextupole;

typedef enum {
  ELETYPE_DRIFT,
  ELETYPE_QUAD,
  ELETYPE_SBEND,
  ELETYPE_MULTIPOLE,
  ELETYPE_SEXTUPOLE,
} EleType;

typedef struct {
  EleType type;
  double R_matrix[BEAM_DOFS*BEAM_DOFS];
  union {
    Drift drift;
    Quad quad;
    Sbend sbend;
    Multipole multipole;
    Sextupole sextupole;
  } as;
} Element;

typedef struct {
  char* key;
  Element value;
} ElementLibrary;

double synch_rad_integral_1(Element *line, size_t periodicity);
double synch_rad_integral_2(Element *line, size_t periodicity);
double synch_rad_integral_3(Element *line, size_t periodicity);
double element_length(Element element);
double bending_radius_of_element(Element element);
Element create_element(char **cursor);
double assigned_double_from_string(char *string);
void element_print(Element element);
bool isvalididchar(char c);
char *populate_element_library(char *cursor);
double calculate_line_length(Element *line);
double calculate_line_angle(Element *line);
void create_line(char *cursor, Element **line);
float e_loss_per_turn(float I2, float gamma0);

#endif // !ELEMENTS_H

