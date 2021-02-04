#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <ctype.h>

#include "utils.h"

/// Check if provided strings are equal
bool str_compare(const char *str_1, const char *str_2) {
  int ret = strcmp(str_1, str_2);
  return (ret == 0) ? true : false;
}

/// If given pointer is null report problem and optionally exit from the program
bool null_report_problem_and_optionally_exit(void *ptr, bool should_exit, const char* name_of_function, const char *msg) {
  if (ptr == NULL) {
    set_text_red();
    printf("\n[Error]:");
    reset_text_color(); 
    printf("(occurred at %s): %s", name_of_function, msg);
    if (should_exit) {
      exit(1);
    }
    return true;
  }
  return false;
}

/// If error occurred report problem and optionally exit from the program
bool error_occurred_report_problem_and_optionally_exit(bool error_occurred, bool should_exit, const char* name_of_function, const char *msg) {
  if (error_occurred) {
    set_text_red();
    printf("\n[Error]: ");
    reset_text_color();
    printf("(occurred at %s): %s", name_of_function, msg);
    if (should_exit) {
      exit(1);
    }
    return true;
  }
  return false;
}

/// Report problem to a user and exit if 'should_exit' parameter set
/// @param should_exit whether error should terminate program
/// @param message message to print to a user
void report_problem_to_user(bool should_exit, const char *message) {
  set_text_red();
  printf("[Error]: ");
  reset_text_color();

  printf("%s\n", message);

  if (should_exit) {
    exit(1);
  }
}

/// Report error when found unknown parameter in configuration file
void report_unknown_parameter(const char* name_of_parameter, uint32_t number_of_line) {
  set_text_red();
  printf("[Error]: ");
  reset_text_color();

  printf("Unknown configuration parameter '%s' in config file at line %u\n", name_of_parameter, number_of_line);
}

/// Report error when found unknown value of a known paramter in configration file
void report_unknown_parameter_value(const char *name_of_parameter, const char* name_of_value, uint32_t number_of_line) {
  set_text_red();
  printf("[Error]: ");
  reset_text_color();

  printf("Unknown value '%s' of parameter '%s' in config file at line %u\n", name_of_value, name_of_parameter, number_of_line);
}

size_t get_size_of_file(FILE *file_stream) {
  // Determine size of the configuration file to allocate enough memory for it
  fseek(file_stream, 0, SEEK_END);
  size_t size_of_file = ftell(file_stream);
  // Return to the start of the configuration file
  fseek(file_stream, 0, SEEK_SET);

  return size_of_file;
}

/// Trim non-alphanumeric symbols at the start and the end
/// @return Newly allocated trimmed line
char* trim_non_alphanumeric_start_end(const char *line) {
  size_t len = strlen(line);

  // Allocate memory for trimmed line
  char *trimmed_line = (char*)allocate_and_zero(len + 1);

  // Add to the new line only alphanumeric characters
  for (size_t orig_i = 0, res_i = 0; orig_i < len; orig_i++) {
    uint8_t ch = line[orig_i];
    if (isalnum(ch) || ch == '=' || ch == '\0' || ch == '_' || ch == '-' || ch == '.' || ch == '/' || ch == '\\') {
      trimmed_line[res_i] = ch;
      res_i++;
    }
  }
  trimmed_line[strlen(trimmed_line)] = '\0';

  // Deallocate original string
  deallocate_and_null((void**)&line);

  return trimmed_line;
}

/// Sets following text color to red
void set_text_red () {
  printf("\033[1;31m");
}

/// Sets following text color to yellow
void set_text_yellow() {
  printf("\033[1;33m");
}

/// Sets following text color to green
void set_text_green() {
  printf("\033[0;32m");
}

/// Resets following text color to default
void reset_text_color () {
  printf("\033[0m");
}

/// Allocates memory and zeros it
/// @param len length of allocation block in bytes
void* allocate_and_zero(size_t len) {
  // Non-zero length
  assert(len > 0 && "Got zero length");

  // Allocate memory and zero it
  void* ptr = calloc(len, sizeof(uint8_t));
  null_report_problem_and_optionally_exit(ptr, true, "allocate_and_zero", "Cannot allocate memory");

  return ptr;
}

/// Deallocates memory and nulls a pointer
/// @param ptr pointer to be deallocated
void deallocate_and_null(void** ptr) {
  null_report_problem_and_optionally_exit(*ptr, true, "deallocate_and_null", "Trying to free null pointer");
  
  free(*ptr);
  *ptr = NULL;
}