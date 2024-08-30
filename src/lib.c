#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <assert.h>

#include "lib.h"

#define PI 3.14159265358979323846

Arena make_arena(void) {
  Arena result = {0};
  result.capacity = ARENA_CAP;
  result.data = calloc(result.capacity, 1);
  if (result.data == NULL) {
    fprintf(stderr, "Memory issues. Dying.\n");
    exit(1);
  }
  return result;
}

void arena_free(Arena *arena) {
  free(arena->data);
  arena->capacity = 0;
  arena->size = 0;
}

void *arena_alloc(Arena *arena, size_t size) {
  if (arena->size + size > arena->capacity) {
    fprintf(stderr, "Memory issues. Dying.\n");
    exit(1);
  }
  void *result = arena->data + arena->size;
  arena->size += size;
  return result;
}

double radians_to_degrees(double radians) {
  return radians * 180.0f / PI;
}

char *join_lines(char* cursor) {
  char *result = cursor;

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

  return result;
}

size_t read_entire_file_to_lines(char *file_path, char **buffer, char ***lines) {
  *buffer = read_entire_file(file_path);
  size_t num_lines = string_to_lines(buffer, lines);
  return num_lines;
}

char *read_entire_file(char *file_path) {
  // Reads an entire file into a char array, and returns a ptr to this. The ptr should be freed by the caller
  FILE *f = fopen(file_path, "r");
  if (f==NULL) {
    fprintf(stderr, "Could not read %s: %s\n", file_path, strerror(errno));
    exit(1);
  }

  fseek(f, 0L, SEEK_END);
  size_t sz = (size_t)ftell(f);
  fseek(f, 0L, SEEK_SET);

  char *contents = calloc(2*sz, sizeof(char));
  if (contents==NULL) {
    fprintf(stderr, "Could not allocate memory. Buy more RAM I guess?\n");
    exit(1);
  }
  fread(contents, 1, sz, f);

  fclose(f);
  
  return contents;
}

size_t string_to_lines(char **string, char ***lines) {
  char *cursor = *string;
  size_t num_lines = count_lines(cursor);

  *lines = calloc(num_lines, sizeof(char*));

  size_t line_ctr = 0;
  while (*cursor) {
    (*lines)[line_ctr] = cursor;
    assert(line_ctr < num_lines);
    advance_to_char(&cursor, '\n');
    *cursor = '\0';
    (cursor)++;
    line_ctr++;
  }

  return num_lines;
}

void advance_to_char(char **string, char c) {
  while (**string != c) (*string)++;
}

size_t count_lines(char *contents) {
  size_t result = 0;
  while (*contents) {
    if (*contents == '\n') result++;
    contents++;
  }
  return result;
}

void advance_to_next_line(char **string) {
  while (**string != '\n') (*string)++;
  (*string)++;
}

