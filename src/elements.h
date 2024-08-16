#ifndef ELEMENTS_H
#define ELEMENTS_H

typedef struct {
  float length;
} Drift;

typedef struct {
  float length;
  float k;
} Quad;

typedef enum {
  ELETYPE_DRIFT,
  ELETYPE_QUAD,
} EleType;

typedef struct {
  EleType type;
  union {
    Drift drift;
    Quad quad;
  } as;
} Element;

#endif // !ELEMENTS_H
#define ELEMENTS_H

