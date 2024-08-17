#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <stdbool.h>
#include <string.h>

#include "lib.h"
#include "elements.h"

#define STB_DS_IMPLEMENTATION
#include "stb_ds.h"

#define ARENA_IMPLEMENTATION
#include "arena.h"

typedef struct {
  char* key;
  Element value;
} hash_map;

static Arena memory_arena = {0};

Element create_element(char **cursor) {
  char type[10] = {0};
  Element result = {0};

  while (**cursor == ' ' | **cursor == ':') (*cursor)++;

  size_t i = 0;
  while (isalpha(**cursor)) {
    type[i++] = toupper(**cursor);
    (*cursor)++;
  }

  if (strcmp(type, "DRIFT") == 0) {
    result.type = ELETYPE_DRIFT;
    while (strncmp(*cursor, "L", 1) != 0) (*cursor)++;
    while (!isdigit(**cursor)) (*cursor)++;
    result.as.drift.length = strtof(*cursor, cursor);
  } else if (strcmp(type, "MARKER") == 0 | strcmp(type, "KICKER") == 0) {
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
    char *k3l_string = strcasestr(*cursor, "k3l");
    if (k3l_string) {
      while (!isdigit(*k3l_string)) k3l_string++;
      result.as.multipole.K3L = strtof(k3l_string, NULL);
    }
    while (**cursor != '\0') (*cursor)++;
    (*cursor)++;
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

int main(void) {
  char *file_path = "./lattices/max4u_lattice.mad8";
  char *buffer = read_entire_file(file_path);
  char *cursor = buffer;

  cursor = join_lines(cursor);

  hash_map *hmap = NULL;

  while (*cursor != '\0') {
    if (*cursor == '!' | *cursor == '&') {
      advance_to_next_line(&cursor);
    } else if (isalpha(*cursor)) {
      // char *element_name = cursor;
      while (isvalididchar(*cursor)) {
        cursor++;
      }
      *cursor = '\0';
      cursor++;

      Element element = create_element(&cursor);

      element_print(element);
    } else {
      cursor++;
    }
  }

  shfree(hmap);
  arena_free(&memory_arena);
  free(buffer);

  printf("\n");

  return 0;
}
