#pragma once

#include <inttypes.h>
#include <stdio.h>

#include "graph.h"

bool str_compare(const char *str_1, const char *str_2);
bool null_report_problem_and_optionally_exit(void *ptr, bool should_exit, const char* name_of_function, const char *msg);
bool error_occurred_report_problem_and_optionally_exit(bool error_occurred, bool should_exit, const char* name_of_function, const char *msg);
void report_unknown_parameter(const char* name_of_parameter, uint32_t number_of_line);
void report_unknown_parameter_value(const char* name_of_parameter, const char *name_of_value, uint32_t number_of_line);
void report_problem_to_user(bool should_exit, const char *message);
void* allocate_and_zero(size_t len);
void deallocate_and_null(void** ptr);
char* trim_non_alphanumeric_start_end(const char *line);
size_t get_size_of_file(FILE *file_stream);
bool is_edge_contains_in(Tuple<Edge> *edges, Edge edge, TypeOfGraph graph_type);
double stop_time(uint64_t start_time);

#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_RESET   "\x1b[0m"

void set_text_red();
void set_text_yellow();
void set_text_green();
void reset_text_color();

#define EPSILON 1e-06