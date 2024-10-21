#ifndef _LIB
#define _LIB

#include <stddef.h>
#include <stdbool.h>

#define ARENA_CAP 1024 * 1024

typedef struct {
  size_t capacity;
  size_t size;
  void *data;
} Arena;

typedef struct {
  char* key;
  double value;
} str2dbl_hashmap;

typedef struct {
  size_t periodicity;
  int harmonic_number;
  double E_0;
  char *file_path;
  char *programname;
  char *twiss_filename;
  bool save_twiss;
} CommandLineArgs;

Arena make_arena(void);
void arena_free(Arena *arena);
void *arena_alloc(Arena *arena, size_t size);

CommandLineArgs get_clargs(int argc, char **argv);
bool isvalididchar(char c);
double radians_to_degrees(double radians);
char *join_lines(char* cursor);
size_t read_entire_file_to_lines(char *file_path, char **buffer, char ***lines);
char *read_entire_file(const char *file_path);
size_t string_to_lines(char **string, char ***lines);
void advance_to_char(char **string, char c);
size_t count_lines(char *contents);
void advance_to_next_line(char **string);

#endif // !_LIB

