#ifndef _LIB
#define _LIB

#include <stddef.h>
#include <stdbool.h>

typedef struct {
  char* key;
  double value;
} str2dbl_hashmap;

typedef struct {
  size_t periodicity;
  double E_0;
  char *file_path;
  char *programname;
  char *twiss_filename;
  bool save_twiss;
  bool error;
} CommandLineArgs;

CommandLineArgs get_clargs(int argc, char **argv);
bool isvalididchar(char c);
double degrees_to_radians(double degrees);
double radians_to_degrees(double radians);
char *join_lines(char* cursor);
size_t read_entire_file_to_lines(char *file_path, char **buffer, char ***lines);
char *read_entire_file(const char *file_path);
size_t string_to_lines(char **string, char ***lines);
void advance_to_char(char **string, char c);
size_t count_lines(char *contents);
void advance_to_next_line(char **string);

#endif // !_LIB

