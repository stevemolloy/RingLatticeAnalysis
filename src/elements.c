#include <stdio.h>

#include "elements.h"

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
      printf("Multipole: K3L = %f\n", element.as.multipole.K3L);
      break;
    case ELETYPE_SEXTUPOLE:
      printf("Sextupole: L = %f, K2 = %f\n", 
             element.as.sextupole.length,
             element.as.sextupole.K2);
      break;
  }
}

