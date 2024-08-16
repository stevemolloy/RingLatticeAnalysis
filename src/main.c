#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include "lib.h"
#include "elements.h"

#define STB_DS_IMPLEMENTATION
#include "stb_ds.h"

#define ARENA_IMPLEMENTATION
#include "arena.h"

typedef struct {
  char* key;
  float value;
} hash_map;

static Arena memory_arena = {0};

int main(void) {
  char *file_path = "./lattices/max4u_lattice.mad8";
  char *buffer = read_entire_file(file_path);
  char *cursor = buffer;

  while (*cursor != '\0') {
    if (*cursor == '&' | *cursor == '!') {
      while (*cursor != '\n') {
        *cursor = ' ';
        cursor++;
      }
      *cursor = ' ';
    } else {
      cursor++;
    }
  }

  cursor = buffer;
  printf("%s\n", cursor);

  hash_map *hmap = NULL;

  while (*cursor != '\0') {
    if (*cursor == '!' | *cursor == '&') {
      advance_to_next_line(&cursor);
    } else if (*cursor == '\n' /*|  *cursor == ' ' */) {
      cursor++;
    } else if (isalpha(*cursor)) {
      char *word = cursor;
      while (isalnum(*cursor)) {
        cursor++;
      }
      *cursor = '\0';
      cursor++;
      shput(hmap, arena_strdup(&memory_arena, word), (float)3.14);
    } else {
      cursor++;
    }
  }

  printf("%s has value %f\n", "bpm", shget(hmap, "bpm"));

  shfree(hmap);
  arena_free(&memory_arena);
  free(buffer);

  printf("\n");

  return 0;
}
