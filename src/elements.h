#ifndef ELEMENTS_H
#define ELEMENTS_H

typedef struct {
  float length;
} Drift;

typedef struct {
  float length;
  float k;
} Quad;

// SBEND,L=0.20424000,  ANGLE = 0.00388598, K1 = -1.33342000  E1= 0.00000000,  E2= 0.00000000 ; 
typedef struct {
  float length;
  float angle;
  float K1;
  float E1;
  float E2;
} Sbend;

typedef enum {
  ELETYPE_DRIFT,
  ELETYPE_QUAD,
  ELETYPE_SBEND,
} EleType;

typedef struct {
  EleType type;
  union {
    Drift drift;
    Quad quad;
    Sbend sbend;
  } as;
} Element;

#endif // !ELEMENTS_H

