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

static Arena memory_arena = {0};

bool isvalididchar(char c) {
  return isalnum(c) | (c == '_');
}

Element* element_list = NULL;
ElementLibrary *element_library = NULL;

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

  printf("The following remains to be parsed:\n%s", cursor);

  char *test_ele_name = "d2_3";
  Element test_ele = shget(element_library, test_ele_name);
  printf("Element %s is ", test_ele_name);
  element_print(test_ele);

  arrfree(element_list);
  shfree(element_library);
  arena_free(&memory_arena);
  free(buffer);

  printf("\n");

  return 0;
}
