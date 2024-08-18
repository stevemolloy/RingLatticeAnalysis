#include <ctype.h>
#include <stdbool.h>

#include "lib.h"
#include "elements.h"

#define STB_DS_IMPLEMENTATION
#include "stb_ds.h"

#define ARENA_IMPLEMENTATION
#include "arena.h"

typedef struct {
  char* key;
  Element value;
} ElementLibrary;

bool isvalididchar(char c) {
  return isalnum(c) | (c == '_');
}

Element* element_list = NULL;
ElementLibrary *element_library = NULL;
Element* line = NULL;

int main(void) {
  char *file_path = "./lattices/max4u_lattice.mad8";
  char *buffer = read_entire_file(file_path);
  char *cursor = buffer;

  cursor = join_lines(cursor);

  while (*cursor != '\0') {
    if (*cursor == '!' | *cursor == '&') {
      advance_to_next_line(&cursor);
    } else if (isalpha(*cursor)) {
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

  advance_to_char(&cursor, '(');
  while (*cursor != ')') {
    while (!isalpha(*cursor) & (*cursor != ')')) {
      cursor++;
    }
    char *temp_cursor = cursor;
    while (isvalididchar(*temp_cursor)) {
      temp_cursor++;
    }
    *temp_cursor = '\0';
    Element ele = shget(element_library, cursor);
    arrput(line, ele);
    cursor = temp_cursor + 1;
  }

  float total_length = 0;
  for (size_t i=0; i<arrlenu(line); i++) {
    element_print(line[i]);
    switch (line[i].type) {
      case ELETYPE_QUAD:
        total_length += line[i].as.quad.length;
        break;
      case ELETYPE_DRIFT:
        total_length += line[i].as.drift.length;
        break;
      case ELETYPE_SBEND:
        total_length += line[i].as.sbend.length;
        break;
      case ELETYPE_MULTIPOLE:
        total_length += line[i].as.multipole.length;
        break;
      case ELETYPE_SEXTUPOLE:
        total_length += line[i].as.sextupole.length;
        break;
    }
  }
  printf("Total length of this line is %f m\n", total_length);

  arrfree(line);
  arrfree(element_list);
  shfree(element_library);
  free(buffer);

  printf("\n");

  return 0;
}
