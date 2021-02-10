#include "config_param_parsing.h"

#include <string.h>
#include <stdio.h>
#include <stdlib.h>

/// Parses line and returns configuration parameter struct
/// @return parsed configuration parameter
ConfigParam parse_configuration_param(const char *orig_line)
{
  // Make a copy of original line 'cause it's gonna be freed by the following function
  size_t original_line_length = strlen(orig_line) + 1;
  char *line_copy = (char*)allocate_and_zero(original_line_length);
  memcpy(line_copy, orig_line, original_line_length);
  // Trim all non-alphanumeric symbols
  const char *line = trim_non_alphanumeric_start_end(line_copy);

  // Length of the input line
  size_t len = strlen(line);
  if (len == 0) { return ConfigParam({ NULL, NULL }); }

  // Create configuration parameter
  ConfigParam param;
  param.init(len);

  // Iterate through line and find parameter name and its value
  bool param_name = true;
  for (size_t i = 0; i < len; i++) {
    // Check when parameter name end and starts value of the parameter  
    if (line[i] == '=') {
      param_name = false;
      continue;
    }

    if (param_name) {
      // Write charactares to parameter name
      param.push_to_name(line[i]);
    } else {
      // Write characters to parameter value
      param.push_to_value(line[i]);
    }
  }

  return param;
}

/// Loads configuration parameters and set configuration flags for the program
/// @param number_of_arguments number of arguments of the current process
/// @param arguments given to the current process
ConfigFlags load_config_file(int number_of_arguments, char **arguments) {
  // Check if path to configuration file is provided
  if (number_of_arguments < 2) {
    printf("Name of configuration file is not provided!\n");
    exit(0);
  }

  // Set configuration flags
  return read_and_set_config_flags(arguments[1]);
}

/// Read configuration file and set configuration flags
ConfigFlags read_and_set_config_flags(const char* path_to_config_file) {
  // Prompt
  printf("Reading configuration file\n");

   // Open configuration file
  FILE *config_file = NULL;
  int error_code = fopen_s(&config_file, path_to_config_file, "r");
  if (error_code != 0) {
    report_problem_to_user(true, "Cannot open configuration file");
  }

  // Get total size of the config file
  int size_of_file = get_size_of_file(config_file);

  // Allocate memory for configuration file
  char *line_buffer = (char*)allocate_and_zero(size_of_file);
  null_report_problem_and_optionally_exit(line_buffer, true, "read_and_set_config_flags", "Cannot allocate memory for configuration file buffer!\n");

  // Configuration parameters
  TypeOfGraph type_of_graph = DIRECTED;
  char *path_to_test_file = NULL;
  bool enable_first_neighborhood = false;
  int attempt_limit = 1; 

  // Number of line to report in case of an error
  uint32_t number_of_line = 1;
  // Read configuration file
  while (fgets(line_buffer, size_of_file, config_file) != NULL) {
    // Parse line
    ConfigParam param = parse_configuration_param(line_buffer);

    if (str_compare(param.config_param_name, "graph")) {
      // If this parameter is type of a graph parameter
      if (str_compare(param.config_param_value, "directed")) {
        type_of_graph = DIRECTED;
        printf("Type of graph: directed\n");
      } else if (str_compare(param.config_param_value, "undirected")) {
        printf("Type of graph: undirected\n");
        type_of_graph = UNDIRECTED;
      } else {
        printf("Type of graph: directed\n");
        report_unknown_parameter_value(param.config_param_name, param.config_param_value, number_of_line);
      }
    } else if (str_compare(param.config_param_name, "test_file")) {
      path_to_test_file = (char*)allocate_and_zero(strlen(param.config_param_value) + 1);
      memcpy(path_to_test_file, param.config_param_value, strlen(param.config_param_value) + 1);
    } else if (str_compare(param.config_param_name, "first_neighborhood")) {
      if (str_compare(param.config_param_value, "enable")) {
        enable_first_neighborhood = true;
      } else if (str_compare(param.config_param_value, "disable")) {
        enable_first_neighborhood = false;
      } else {
        printf("First neighborhood: disabled\n");
        report_unknown_parameter_value(param.config_param_name, param.config_param_value, number_of_line);
      }
    } else if (str_compare(param.config_param_name, "attempt_limit")) {
      attempt_limit = std::atoi(param.config_param_value);
      printf("Attempt limit: %d\n", attempt_limit);
    } else {
      // Unkown type of parameter
      report_unknown_parameter(param.config_param_name, number_of_line);
      exit(0);
    }

    // Free memory
    param.free(); 

    number_of_line++;   
  }

  // Close configuration file
  fclose(config_file);
  // Deallocate memory for the configuration file buffer
  deallocate_and_null((void**)&line_buffer);

  return ConfigFlags({ 
    .graph_type = type_of_graph,
    .path_to_test_file = path_to_test_file,
    .first_neighborhood_enabled = enable_first_neighborhood,
    .attempt_limit = attempt_limit,
  });
}