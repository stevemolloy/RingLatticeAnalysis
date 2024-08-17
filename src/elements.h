#ifndef ELEMENTS_H
#define ELEMENTS_H

typedef struct {
  float length;
} Drift;

typedef struct {
  float length;
  float K1;
} Quad;

// SBEND,L=0.20424000,  ANGLE = 0.00388598, K1 = -1.33342000  E1= 0.00000000,  E2= 0.00000000 ; 
typedef struct {
  float length;
  float angle;
  float K1;
  float E1;
  float E2;
} Sbend;

// Multipole, K3L= -38933.64000000; 
typedef struct {
  float K3L;
} Multipole;

//Sextupole,  L= 0.10000000, K2= -272.67400000 ;
typedef struct {
  float length;
  float K2;
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
  union {
    Drift drift;
    Quad quad;
    Sbend sbend;
    Multipole multipole;
    Sextupole sextupole;
  } as;
} Element;

Element create_element(char **cursor);
void element_print(Element element);

#endif // !ELEMENTS_H

