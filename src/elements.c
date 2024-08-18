#define _GNU_SOURCE
#include <string.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>

#include "elements.h"

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
    while (strncmp(*cursor, "L", 1) != 0) (*cursor)++;
    while (!isdigit(**cursor)) (*cursor)++;
    result.as.drift.length = strtof(*cursor, cursor);
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
    if (l_string) {
      while (!isdigit(*l_string)) l_string++;
      result.as.sbend.length = strtof(l_string, NULL);
    }
    if (angle_string) {
      while (!isdigit(*angle_string)) angle_string++;
      result.as.sbend.angle = strtof(angle_string, NULL);
    }
    if (k1_string) {
      while (!isdigit(*k1_string)) k1_string++;
      result.as.sbend.angle = strtof(k1_string, NULL);
    }
    if (e1_string) {
      while (!isdigit(*e1_string)) e1_string++;
      result.as.sbend.angle = strtof(e1_string, NULL);
    }
    if (e2_string) {
      while (!isdigit(*e2_string)) e2_string++;
      result.as.sbend.angle = strtof(e2_string, NULL);
    }
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
    if (l_string) {
      while (!isdigit(*l_string)) l_string++;
      result.as.quad.length = strtof(l_string, NULL);
    }
    if (k1_string) {
      while (!isdigit(*k1_string)) k1_string++;
      result.as.quad.K1 = strtof(k1_string, NULL);
    }
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
    if (l_string) {
      while (!isdigit(*l_string)) l_string++;
      result.as.sextupole.length = strtof(l_string, NULL);
    }
    if (k2_string) {
      while (!isdigit(*k2_string)) k2_string++;
      result.as.sextupole.K2 = strtof(k2_string, NULL);
    }
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
    if (k3l_string) {
      while (!isdigit(*k3l_string)) k3l_string++;
      result.as.multipole.K3L = strtof(k3l_string, NULL);
    }
    if (l_string) {
      while (!isdigit(*l_string)) l_string++;
      result.as.multipole.length = strtof(l_string, NULL);
    }
    while (**cursor != '\0') (*cursor)++;
    (*cursor)++;
  } else if (strcmp(type, "LINE") == 0) {
    *cursor -= 4; // Safe since we just extracted "LINE" from this, so this puts the cursor back there.
    printf("Found a LINE, so returning to the main function\n");
    return result;
  } else {
    fprintf(stderr, "ERROR: Unable to parse elements of type %s\n", type);
    printf("%s\n", *cursor);
    exit(1);
  }

  return result;
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

