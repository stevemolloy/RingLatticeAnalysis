#ifndef _LIB
#define _LIB

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
#define STRERROR strerror_s
#else
#define STRERROR strerror
#endif

#include <stddef.h>

#define ARENA_CAP 1024 * 1024

typedef struct {
  size_t capacity;
  size_t size;
  void *data;
} Arena;

Arena make_arena(void);
void arena_free(Arena *arena);
void *arena_alloc(Arena *arena, size_t size);

double radians_to_degrees(double radians);
char *join_lines(char* cursor);
size_t read_entire_file_to_lines(char *file_path, char **buffer, char ***lines);
char *read_entire_file(const char *file_path);
size_t string_to_lines(char **string, char ***lines);
void advance_to_char(char **string, char c);
size_t count_lines(char *contents);
void advance_to_next_line(char **string);

#endif // !_LIB

