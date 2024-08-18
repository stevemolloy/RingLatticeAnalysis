#define _GNU_SOURCE
#include <string.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>

#include "stb_ds.h"

#include "elements.h"

extern Element* element_list;
extern ElementLibrary *element_library;

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
    result.as.drift.length = assigned_float_from_string(l_string);
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
    result.as.sbend.length = assigned_float_from_string(l_string);
    result.as.sbend.angle = assigned_float_from_string(angle_string);
    result.as.sbend.K1 = assigned_float_from_string(k1_string);
    result.as.sbend.E1 = assigned_float_from_string(e1_string);
    result.as.sbend.E2 = assigned_float_from_string(e2_string);
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
    result.as.quad.length = assigned_float_from_string(l_string);
    result.as.quad.K1 = assigned_float_from_string(k1_string);
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
    result.as.sextupole.length = assigned_float_from_string(l_string);
    result.as.sextupole.K2 = assigned_float_from_string(k2_string);
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
    result.as.multipole.length = assigned_float_from_string(l_string);
    result.as.multipole.K3L = assigned_float_from_string(k3l_string);
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

float assigned_float_from_string(char *string) {
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

