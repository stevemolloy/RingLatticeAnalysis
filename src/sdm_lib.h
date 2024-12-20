#ifndef _SDM_LIB_H
#define _SDM_LIB_H

/* This library provides the following:
 *
 * char *sdm_shift_args(int *argc, char ***argv);      Peel arguments off the **argv array typically provided to main, decrementing argc appropriately.
 * char *sdm_read_entire_file(const char *file_path);  Read the contents of a file into a character array. This character array is malloc'ed and so should be freed by the user.
 * SDM_FREE_AND_NULL(ptr)                              Free the memory pointed to by ptr, and then set ptr to NULL.
 * #define SDM_FREE SDM_FREE_AND_NULL
 * #define SDM_MALLOC malloc
 * 
 * # DYNAMIC ARRAYS
 * ================
 * SDM_ENSURE_ARRAY_CAP(da, cap)       Ensure that the dynamic array, da, has a capacity equal-to cap. Realloc is used if needed.
 * SDM_ENSURE_ARRAY_MIN_CAP(da, cap)   Ensure that the dynamic array, da, has a capacity equal-to or greater-than cap. Realloc is used if needed.
 * DEFAULT_CAPACITY 128                Default capacity in bytes to be used when the capacity has not been set by the user.
 * SDM_ARRAY_PUSH(da, item)            Push the value of item to the dynamic array, da, reallocing if needed.
 * SDM_ARRAY_SWAP(da, ind1, ind2)      Swap the elements at the marked indices (if length==capacity this will extend the array)
 * SDM_ARRAY_FREE(da)                  Free the memory in the dynamic array, da, and zero things
 * SDM_ARRAY_RESET(da)                 Reset the length of the dynamic array, da, to zero, effectively emptying it. No memory is freed by this.
 * 
 * # STRING VIEWS
 * ==============
 * sdm_string_view sdm_cstr_as_sv(char *cstr);                                 Return a string view containing the provided c-string
 * sdm_string_view sdm_sized_str_as_sv(char *cstr, size_t length);             Return a string view containing the c-string with the provided size.
 * sdm_string_view sdm_sv_pop_by_delim(sdm_string_view *SV, const char delim); Pop off the start of the provided string view up to the given delimiter.
 * void sdm_sv_trim(sdm_string_view *SV);                                      Trim any whitespace from the start of the string view
 * SDM_SV_F "%.*s"                                                             A printf helper.
 * SDM_SV_Vals(S) (int)(S).length, (S).data                                    A printf helper.
 * 
 * # MEMORY ARENA
 * ==============
 * #define SDM_ARENA_DEFAULT_CAP 256 * 1024*1024              Default capacity of the memory arena when not supplied by the user
 * void sdm_arena_init(sdm_arena_t *arena, size_t capacity);  Initialise a memory arena with a certain capacity and malloc the required space.
 * void *sdm_arena_alloc(sdm_arena_t *arena, size_t size);    Allocate a region of size bytes in the given arena, and return a pointer to the start of this region.
 * void sdm_arena_free(sdm_arena_t *arena);                   Deallocate all memory in the arena, and zero everything
 */

#include <stdbool.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#define SDM_FREE_AND_NULL(ptr) \
do {                           \
  if ((ptr)) free((ptr));      \
  (ptr) = NULL;                \
} while (0)

#define SDM_FREE SDM_FREE_AND_NULL
#define SDM_MALLOC malloc

#define hash jenkins_one_at_a_time_hash

#define SDM_ENSURE_ARRAY_CAP(da, cap) do {                     \
    (da).capacity = cap;                                       \
    (da).data = realloc((da).data,                             \
        (da).capacity * sizeof((da).data[0]));                 \
    if ((da).data == NULL) {                                   \
      fprintf(stderr, "ERR: Couldn't alloc memory.\n");        \
      exit(1);                                                 \
    }                                                          \
  } while (0)

#define SDM_ENSURE_ARRAY_MIN_CAP(da, cap) do {                 \
    if ((da).capacity < cap) {                                 \
      (da).capacity = cap;                                     \
      (da).data = realloc((da).data,                           \
          (da).capacity * sizeof((da).data[0]));               \
      if ((da).data == NULL) {                                 \
        fprintf(stderr, "ERR: Couldn't alloc memory. \n");     \
        exit(1);                                               \
      }                                                        \
    }                                                          \
  } while (0)

#define SDM_ARRAY_POP(da) (da).data[--((da).length)]

#define DEFAULT_CAPACITY 128

#define SDM_ARRAY_PUSH(da, item) do {                             \
    if (((da).capacity == 0) || ((da).data == NULL)) {            \
      (da).capacity = DEFAULT_CAPACITY;                           \
      (da).data = malloc((da).capacity * sizeof((da).data[0]));   \
      memset((da).data, 0, (da).capacity * sizeof((da).data[0])); \
      if ((da).data == NULL) {                                    \
        fprintf(stderr, "ERR: Couldn't alloc memory.\n");         \
        exit(1);                                                  \
      }                                                           \
    }                                                             \
    while ((da).length >= (da).capacity) {                        \
      (da).capacity *= 2;                                         \
      (da).data = realloc((da).data,                              \
           (da).capacity * sizeof((da).data[0]));                 \
      if ((da).data == NULL) {                                    \
        fprintf(stderr, "ERR: Couldn't alloc memory.\n");         \
        exit(1);                                                  \
      }                                                           \
    }                                                             \
    (da).data[(da).length++] = item;                              \
  } while (0);

#define SDM_ARRAY_SWAP(da, ind1, ind2)                                   \
do {                                                                     \
  assert((int)ind1 < (int)da.length && "First index is out of bounds");  \
  assert((int)ind2 < (int)da.length && "Second index is out of bounds"); \
  SDM_ARRAY_PUSH(da, da.data[ind1]);                                     \
  da.data[ind1] = da.data[ind2];                                         \
  da.data[ind2] = da.data[da.length - 1];                                \
  da.length--;                                                           \
} while (0)

#define SDM_ARRAY_FREE(da) do {                                \
    SDM_FREE((da).data);                                           \
    (da).length = 0;                                           \
    (da).capacity = 0;                                         \
  } while (0);

#define SDM_ARRAY_RESET(da) do { (da).length = 0; } while (0)

typedef struct {
  size_t length;
  char *data;
} sdm_string_view;

#define SDM_SV_F "%.*s"
#define SDM_SV_Vals(S) (int)(S).length, (S).data

char *sdm_shift_args(int *argc, char ***argv);

char *sdm_read_entire_file(const char *file_path);

sdm_string_view sdm_cstr_as_sv(char *cstr);
char *sdm_sv_to_cstr(sdm_string_view sv);
sdm_string_view sdm_sized_str_as_sv(char *cstr, size_t length);
sdm_string_view sdm_pop_by_length(sdm_string_view *SV, size_t len);
sdm_string_view sdm_sv_pop_by_delim(sdm_string_view *SV, const char delim);
sdm_string_view sdm_sv_pop_by_whitespace(sdm_string_view *SV);
int sdm_sv_pop_integer(sdm_string_view *SV);
char sdm_sv_pop_one_char(sdm_string_view *SV);
void sdm_sv_trim(sdm_string_view *SV);

typedef struct {
  char key[32];
  double value;
  bool occupied;
} Dbl_KVPair;

typedef struct {
  Dbl_KVPair *data;
  size_t length;
  size_t capacity;
} DblArray;

#define DEFAULT_HM_CAP 256

#define SET_HM_CAPACITY(hm, cap)                                              \
  do {                                                                        \
    (hm)->capacity = (cap);                                                   \
    if ((hm)->data == NULL) {                                                 \
      (hm)->data = malloc(hm->capacity * sizeof(hm->data[0]));                \
    }                                                                         \
    else {                                                                    \
      (hm)->data = realloc((hm)->data, (hm)->capacity * sizeof(hm->data[0])); \
    }                                                                         \
    if ((hm)->data == NULL) {                                                 \
      fprintf(stderr, "ERR: Can't alloc.\n");                                 \
      exit(1);                                                                \
    }                                                                         \
    memset((hm)->data, 0, (hm)->capacity * sizeof(hm->data[0]));              \
  } while (0)

#define FREE_HASHMAP(hm) \
  do {                   \
    (hm).capacity = 0;   \
    SDM_FREE((hm).data);     \
  } while (0);

#define GET_HASHMAP_INDEX(hm, value_of_key, index_addr)                      \
  do {                                                                       \
    *(index_addr) = -1;                                                      \
    uint32_t location = get_hashmap_location((value_of_key), (hm).capacity); \
    for (size_t offset=0; offset<(hm).capacity; offset++) {                  \
      int new_location = (location + offset) % (hm).capacity;                \
      if ((hm).data[new_location].occupied &&                                \
        (strcmp((hm).data[new_location].key, (value_of_key))==0)) {          \
        *(index_addr) = new_location;                                        \
        break;                                                               \
      }                                                                      \
    }                                                                        \
  } while (0)

#define HM_VAL_AT(hm, index) \
  (hm).data[index].value

#define PUSH_TO_HASHMAP(hm1, value_of_key, value_of_value)                                   \
  do {                                                                                       \
    uint32_t location;                                                                       \
    if ((hm1)->capacity == 0) {                                                              \
      SET_HM_CAPACITY((hm1), DEFAULT_HM_CAP);                                                \
    }                                                                                        \
    if ((hm1)->length >= (hm1)->capacity) {                                                  \
      fprintf(stderr, "Hashmap has insufficient capacity.\n");                               \
      exit(1);                                                                               \
    }                                                                                        \
    location = get_hashmap_location((value_of_key), (hm1)->capacity);                        \
                                                                                             \
    if (!(hm1)->data[location].occupied ||                                                   \
      ((hm1)->data[location].occupied && strcmp((value_of_key), (hm1)->data[location].key)==0)) {  \
      strcpy((hm1)->data[location].key, (value_of_key));                                     \
      (hm1)->data[location].value = (value_of_value);                                        \
      (hm1)->data[location].occupied = true;                                                 \
      (hm1)->length++;                                                                       \
    } else {                                                                                 \
      for (size_t i=1; i<(hm1)->capacity; i++) {                                             \
        size_t new_location = (location + i) % (hm1)->capacity;                              \
        if (!(hm1)->data[new_location].occupied ||                                           \
          ((hm1)->data[new_location].occupied && strcmp((value_of_key), (hm1)->data[new_location].key)==0)) { \
          strcpy((hm1)->data[new_location].key, (value_of_key));                             \
          (hm1)->data[new_location].value = (value_of_value);                                \
          (hm1)->data[new_location].occupied = true;                                         \
          (hm1)->length++;                                                                   \
          break;                                                                             \
        }                                                                                    \
      }                                                                                      \
    }                                                                                        \
  } while (0)

void push_to_dblarray(DblArray *hm, char *key, double value);
uint32_t get_hashmap_location(const char* key, size_t capacity);
uint32_t jenkins_one_at_a_time_hash(const uint8_t* key, size_t length);

#define SDM_ARENA_DEFAULT_CAP 256 * 1024*1024

typedef struct sdm_arena_t sdm_arena_t;

struct sdm_arena_t {
  size_t length;
  size_t capacity;
  void *start;
  sdm_arena_t *next;
};

void sdm_arena_init(sdm_arena_t *arena, size_t capacity);
void *sdm_arena_alloc(sdm_arena_t *arena, size_t size);
void sdm_arena_free(sdm_arena_t *arena);

#endif /* ifndef _SDM_LIB_H */

#ifdef SDM_LIB_IMPLEMENTATION

#include <ctype.h>
#include <errno.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

char *sdm_read_entire_file(const char *file_path) {
  // Reads an entire file into a char array, and returns a ptr to this. The ptr should be freed by the caller
  FILE *f = fopen(file_path, "r");
  if (f==NULL) {
    fprintf(stderr, "Could not read %s: %s\n", file_path, strerror(errno));
    exit(1);
  }

  fseek(f, 0L, SEEK_END);
  int sz = ftell(f);
  fseek(f, 0L, SEEK_SET);

  char *contents = calloc(sz + 1, sizeof(char));
  if (contents==NULL) {
    fprintf(stderr, "Could not allocate memory. Buy more RAM I guess?\n");
    exit(1);
  }
  fread(contents, 1, sz, f);

  fclose(f);
  
  return contents;
}

sdm_string_view sdm_cstr_as_sv(char *cstr) {
  return (sdm_string_view){
    .data = cstr,
    .length = strlen(cstr)
  };
}

char *sdm_sv_to_cstr(sdm_string_view sv) {
  char *ret = malloc(sv.length + 1);
  memset(ret, 0, sv.length + 1);
  memcpy(ret, sv.data, sv.length);
  return ret;
}

sdm_string_view sdm_sized_str_as_sv(char *cstr, size_t length) {
  return (sdm_string_view) {
    .data = cstr,
    .length = length
  };
}

sdm_string_view sdm_pop_by_length(sdm_string_view *SV, size_t len) {
  sdm_string_view retval = (sdm_string_view) {
    .data = SV->data,
    .length = len,
  };
  SV->data += len;
  SV->length -= len;
  return retval;
}

sdm_string_view sdm_sv_pop_by_whitespace(sdm_string_view *SV) {
  sdm_string_view ret = {0};
  ret.data = SV->data;

  while (!isspace(*SV->data) && (SV->length>0)) {
    SV->data++;
    SV->length--;
    ret.length++;
  }
  SV->data++;
  SV->length--;

  return ret;
}

sdm_string_view sdm_sv_pop_by_delim(sdm_string_view *SV, const char delim) {
  sdm_string_view ret = {0};
  ret.data = SV->data;

  while ((*SV->data != delim) && (SV->length>0)) {
    SV->data++;
    SV->length--;
    ret.length++;
  }
  SV->data++;
  SV->length--;

  return ret;
}

int sdm_sv_pop_integer(sdm_string_view *SV) {
  char *new_pos;
  int retval = strtol(SV->data, &new_pos, 0);
  size_t diff = new_pos - SV->data;
  SV->data = new_pos;
  SV->length -= diff;
  return retval;
}

char sdm_sv_pop_one_char(sdm_string_view *SV) {
  char retval = SV->data[0];
  SV->data++;
  SV->length--;
  return retval;
}

void sdm_sv_trim(sdm_string_view *SV) {
  while (isspace(*SV->data) && SV->length>0) {
    SV->data++;
    SV->length--;
  }
}

char *sdm_shift_args(int *argc, char ***argv) {
  if (*argc <= 0) return NULL;
  (*argc)--;
  char **ret = *argv;
  (*argv)++;
  return *ret;
}

uint32_t get_hashmap_location(const char* key, size_t capacity) {
  uint32_t key_hash = hash((uint8_t*)key, strlen(key));
  return key_hash % capacity;
}

void resize_dblarray(DblArray *hm) {
  DblArray resized_array = {0};
  SET_HM_CAPACITY((&resized_array), hm->capacity * 2);
  for (size_t i=0; i<hm->capacity; i++) {
    if (!hm->data[i].occupied) continue;
    push_to_dblarray(&resized_array, hm->data[i].key, hm->data[i].value);
  }
  hm->capacity = resized_array.capacity;
  SDM_FREE(hm->data);
  hm->data = resized_array.data;
}

void push_to_dblarray(DblArray *hm, char *key, double value) {
  if (hm->length >= hm->capacity) resize_dblarray(hm);
  PUSH_TO_HASHMAP(hm, key, value);
}

// https://en.wikipedia.org/wiki/Jenkins_hash_function
uint32_t jenkins_one_at_a_time_hash(const uint8_t* key, size_t length) {
  size_t i = 0;
  uint32_t hash = 0;
  while (i != length) {
    hash += key[i++];
    hash += hash << 10;
    hash ^= hash >> 6;
  }
  hash += hash << 3;
  hash ^= hash >> 11;
  hash += hash << 15;
  return hash;
}

void sdm_arena_init(sdm_arena_t *arena, size_t capacity) {
  arena->start = malloc(capacity);
  if (arena->start == NULL) {
    fprintf(stderr, "Memory problem. Aborting.\n");
    exit(1);
  }
  memset(arena->start, 0, capacity);
  arena->next = malloc(sizeof(*arena->next));
  if (arena->next == NULL) {
    fprintf(stderr, "Memory problem. Aborting.\n");
    exit(1);
  }
  memset(arena->next, 0, sizeof(*arena->next));
  arena->capacity = capacity;
  arena->length = 0;
}

void *sdm_arena_alloc(sdm_arena_t *arena, size_t size) {
  if (arena->start == NULL) {
    size_t capacity = (arena->capacity > 0) ? arena->capacity : SDM_ARENA_DEFAULT_CAP;
    while (capacity < size) {
      capacity *= 2;
    }
    sdm_arena_init(arena, capacity);
  }

  if (arena->capacity - arena->length < size) {
    if (arena->next->capacity == 0) arena->next->capacity = arena->capacity;
    return sdm_arena_alloc(arena->next, size);
  }

  void *return_val = (char*)arena->start + arena->length;

  if (size > 0)
    arena->length += size;
  else
    arena->length += 1;

  return return_val;
}

void sdm_arena_free(sdm_arena_t *arena) {
  if (arena->next) {
    sdm_arena_free(arena->next);
  }

  SDM_FREE_AND_NULL(arena->start);

  arena->length = 0;
  arena->capacity = 0;

  SDM_FREE_AND_NULL(arena->next);
}

#endif // SDM_LIB_IMPLEMENTATION

