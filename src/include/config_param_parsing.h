#pragma once

#include <assert.h>
#include <inttypes.h>
#include "utils.h"
#include "graph.h"

struct ConfigFlags {
  TypeOfGraph graph_type;
  const char *path_to_test_file;
  bool first_neighborhood_enabled;
  int attempt_limit = 1;
  bool second_neighborhood_enabled;
};

/// Struct stores information about configuration parameter
struct ConfigParam
{
  char *config_param_name = NULL;
  size_t size_of_param_name = 0;

  char *config_param_value = NULL;
  size_t size_of_param_value = 0;

  size_t line_len = 0;

  // Add character to the name field
  void push_to_name(char ch) {
    config_param_name[size_of_param_name] = ch; 
    size_of_param_name++;
    assert(size_of_param_name < line_len);
  }

  // Add character to the value field
  void push_to_value(char ch) {
    config_param_value[size_of_param_value] = ch;
    size_of_param_value++;
    assert(size_of_param_value < line_len);
  }

  // Initialize and allocate memory for name and value
  void init(const size_t len) {
    // Allocate memory for configuration parameter name
    config_param_name = (char*)allocate_and_zero(len);
    null_report_problem_and_optionally_exit((void *)config_param_name, true, "ConfigParam::init", "Cannot allocate memory for configuration parameter name");

    // Allocate memory for configuration paramter value
    config_param_value = (char*)allocate_and_zero(len);
    null_report_problem_and_optionally_exit((void *)config_param_value, true, "ConfigParam::init", "Cannot allocate memory for configuration parameter value");

    line_len = len;
  }

  // Deallocate memory 
  void free() {
    deallocate_and_null((void**)&config_param_name);
    deallocate_and_null((void**)&config_param_value);
  }
};

/// Functions 
ConfigParam parse_configuration_param(const char *line);
ConfigFlags read_and_set_config_flags(const char *path_to_config_file);
void report_problem_to_user(bool should_exit, const char *message);
ConfigFlags load_config_file(int number_of_arguments, char **arguments);