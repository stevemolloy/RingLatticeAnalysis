#include <stddef.h>

size_t read_entire_file_to_lines(char *file_path, char **buffer, char ***lines);
char *read_entire_file(char *file_path);
size_t string_to_lines(char **string, char ***lines);
void advance_to_char(char **string, char c);
size_t count_lines(char *contents);
void advance_to_next_line(char **string);

