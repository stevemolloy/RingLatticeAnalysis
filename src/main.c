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

typedef struct {
  size_t len;
  size_t cap;
  char *data;
} SizedString;

SizedString new_sized_string(void) {
  SizedString result = {0};
  result.cap = 256;
  result.data = malloc(result.cap * sizeof(char));
  if (result.data == NULL) result.cap = 0;
  return result;
}

void free_sized_string(SizedString s_string) {
  s_string.cap = 0;
  s_string.len = 0;
  free(s_string.data);
}

int append_to_sized_string(SizedString *ca, char new_element) {
  while (ca->len >= ca->cap) {
    ca->cap *= 2;
    ca->data = realloc(ca->data, ca->cap);
  }
  if (ca->data == NULL) return -1;
  ca->data[ca->len] = new_element;
  ca->len++;
  return 1;
}

static Arena memory_arena = {0};

int main(void) {
  char *file_path = "./lattices/max4u_lattice.mad8";
  char *buffer = read_entire_file(file_path);

  char *cursor = buffer;

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
