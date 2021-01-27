#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string>
#include <sstream>

#include "utils.h"
#include "config_param_parsing.h"
#include "win_console_color.h"
#include "graph.h"
#include "cycle_extraction.h"

#include "scip/scip.h"
#include "scip/scipdefplugins.h"

bool are_edges_the_same_undirectional_graphs(Edge edge1, Edge edge2) {
  return ((edge1.start.number == edge2.start.number) && (edge1.end.number == edge2.end.number)) || 
         ((edge1.start.number == edge2.end.number) && (edge1.end.number == edge2.start.number));
}

/// Compare two edges for equality
/// @param edge1 first edge
/// @param edge2 second edge
/// @return bool whether edges are equal
bool are_edges_the_same_directional_graphs(Edge edge1, Edge edge2) {
  return (edge1.start.number == edge2.start.number) && (edge1.end.number == edge2.end.number);
}

/// Add first constraint to the problem 
/// SUM_(e∈E) X_e = |V|
/// @param scip problem variable
/// @param first Graph of the first cycle
/// @param second Graph of the second cycle
/// @return SCIP_RETCODE
SCIP_RETCODE add_first_constraint(SCIP *scip, Graph first, Graph second) {
  SCIP_CONS *constraint = NULL;
  SCIP_VAR **constraint_variables = NULL;
  
  for (size_t i = 0; i < first.number_of_edges(); i++) {
    sb_push(constraint_variables, first.edges[i].var);
  }

  for (size_t i = 0; i < second.number_of_edges(); i++) {
    sb_push(constraint_variables, second.edges[i].var);
  }

  assert((size_t)sb_count(constraint_variables) == first.number_of_edges() + second.number_of_edges());
  
  SCIP_Real *first_constraint_coefficients = NULL;
  for (int i = 0; i < sb_count(constraint_variables); i++) {
    sb_push(first_constraint_coefficients, 1.0);
  }

  SCIP_CALL(SCIPcreateConsBasicLinear(scip, &constraint, "First constraint", sb_count(constraint_variables), 
      constraint_variables, first_constraint_coefficients, first.number_of_verticies(), first.number_of_verticies()));
  SCIP_CALL(SCIPaddCons(scip, constraint));

  // TODO: remove this print 
  SCIP_CALL(SCIPprintCons(scip, constraint, NULL));

  SCIP_CALL(SCIPreleaseCons(scip, &constraint));

  return SCIP_OKAY;
}

/// Add second constraint to the problem
SCIP_RETCODE add_second_constraint(SCIP *scip, Graph first, Graph second, TypeOfGraph graph_type) {
  if (graph_type == DIRECTED) {
    // Branch for directed graph
    for (size_t i = 0; i < first.number_of_edges(); i++) {
      size_t first_cycle_edge_index = i;
      size_t first_cycle_next_edge_index = ((i + 1) == first.number_of_edges()) ? 0 : i + 1;

      Edge first_cycle_edge = first.edges[first_cycle_edge_index];
      Edge first_cycle_next_edge = first.edges[first_cycle_next_edge_index];
      for (size_t j = 0; j < second.number_of_edges(); j++) {
        size_t second_cycle_edge_index = j;
        size_t second_cycle_next_edge_index = ((j + 1) == second.number_of_edges()) ? 0 : j + 1;

        Edge second_cycle_edge = second.edges[second_cycle_edge_index];
        Edge second_cycle_next_edge = second.edges[second_cycle_next_edge_index];
        if ((first_cycle_edge.end.number == second_cycle_edge.end.number) && (first_cycle_next_edge.start.number == second_cycle_next_edge.start.number)) {
          SCIP_CONS* constraint_ingoing_edges;
          SCIP_CONS* constraint_outgoing_edges;

          SCIP_VAR *variables_ingoing_edges[2] = { };
          SCIP_VAR *variables_outgoing_edges[2] = { };

          Edge first_cycle_ingoing = first_cycle_edge;
          Edge first_cycle_outgoing = first_cycle_next_edge;
          Edge second_cycle_ingoing = second_cycle_edge;
          Edge second_cycle_outgoing = second_cycle_next_edge;

          variables_ingoing_edges[0] = first_cycle_ingoing.var;
          variables_ingoing_edges[1] = second_cycle_ingoing.var;

          variables_outgoing_edges[0] = first_cycle_outgoing.var;
          variables_outgoing_edges[1] = second_cycle_outgoing.var;
        
          SCIP_Real coefficients[2] = { 1.0, 1.0 };

          // Ingoing edges
          SCIP_CALL(SCIPcreateConsBasicLinear(scip, &constraint_ingoing_edges, "Ingoing edges", 2, variables_ingoing_edges, coefficients, 1.0, 1.0));
          SCIP_CALL(SCIPaddCons(scip, constraint_ingoing_edges));

          // TODO: remove this print
          SCIP_CALL(SCIPprintCons(scip, constraint_ingoing_edges, NULL));
          printf("\n");

          SCIP_CALL(SCIPreleaseCons(scip, &constraint_ingoing_edges));

          // Outgoing edges
          SCIP_CALL(SCIPcreateConsBasicLinear(scip, &constraint_outgoing_edges, "Ingoing edges", 2, variables_outgoing_edges, coefficients, 1.0, 1.0));
          SCIP_CALL(SCIPaddCons(scip, constraint_outgoing_edges));
          
          // TODO: remove this print
          SCIP_CALL(SCIPprintCons(scip, constraint_outgoing_edges, NULL));

          SCIP_CALL(SCIPreleaseCons(scip, &constraint_outgoing_edges));

          // Go to the next vertex
          break;
        }
      }
    }
  } else {
    // Branch for undirected graphs
    // 4 is number of incident edges
    SCIP_VAR *constraint_variables[4] = { };

    SCIP_CONS *constraint;

    for (size_t i = 0; i < first.number_of_edges(); i++) {
      size_t first_cycle_edge_index = i;
      size_t first_cycle_next_edge_index = ((i + 1) == first.number_of_edges()) ? 0 : i + 1;

      Edge first_cycle_edge = first.edges[first_cycle_edge_index];
      Edge first_cycle_next_edge = first.edges[first_cycle_next_edge_index];
      for (size_t j = 0; j < second.number_of_edges(); j++) {
        size_t second_cycle_edge_index = j;
        size_t second_cycle_next_edge_index = ((j + 1) == second.number_of_edges()) ? 0 : j + 1;

        Edge second_cycle_edge = second.edges[second_cycle_edge_index];
        Edge second_cycle_next_edge = second.edges[second_cycle_next_edge_index];
        if ((first_cycle_edge.end.number == second_cycle_edge.end.number) && (first_cycle_next_edge.start.number == second_cycle_next_edge.start.number)) {
          constraint_variables[0] = first_cycle_edge.var;
          constraint_variables[1] = first_cycle_next_edge.var;
          constraint_variables[2] = second_cycle_edge.var;
          constraint_variables[3] = second_cycle_next_edge.var;

          SCIP_Real coefficients[] =  { 1.0, 1.0, 1.0, 1.0 };
          SCIP_CALL(SCIPcreateConsBasicLinear(scip, &constraint, "Second constraint", 4, constraint_variables, coefficients, 2.0, 2.0));
          SCIP_CALL(SCIPaddCons(scip, constraint));

          // TODO: remove this print
          SCIP_CALL(SCIPprintCons(scip, constraint, NULL));

          SCIP_CALL(SCIPreleaseCons(scip, &constraint));

          // Move to the next vertex
          break;
        }
      }
    }
  }
  
  return SCIP_OKAY;
}

/// Initialize SCIP
/// @param scip SCIP varable holding problem formulation
/// @return SCIP_RETCODE
SCIP_RETCODE init_scip(SCIP **scip) {
  SCIP_CALL(SCIPcreate(scip));
  SCIP_CALL(SCIPincludeDefaultPlugins(*scip));
  SCIP_CALL(SCIPsetEmphasis(*scip, SCIP_PARAMEMPHASIS_DEFAULT, true));

  SCIP_CALL(SCIPcreateProbBasic(*scip, "ILP"));

  return SCIP_OKAY;
}

int main(int argc, char **argv) {
  setupConsole();

  // Load configuration file into the program
  ConfigFlags flags = load_config_file(argc, argv);

  const char *path_to_the_file = "test_cases\\test_6_undir.txt";

  FILE *file = NULL;
  int error_code = fopen_s(&file, path_to_the_file, "r");
  if (error_code != 0) {
    report_problem_to_user(true, "Couldn't open the file with cycles");
  }

  // Load cycles to the program
  size_t total_file_size = get_size_of_file(file);
  BundleOfGraphs b_graphs = read_cycles_from_stream(file, total_file_size);

  // Set graphs to the corresponding cycles
  Graph first = b_graphs.first;
  Graph second = b_graphs.second;

  SCIP *scip = NULL;

  // Initialize SCIP
  init_scip(&scip);

  // Add variables to the problem from cycles
  for (size_t i = 0; i < first.number_of_verticies(); i++) {
    std::string name_of_variable_first_cycle = (std::string("Fx") + std::to_string(first.edges[i].start.number) + std::to_string(first.edges[i].end.number));
    SCIP_CALL(SCIPcreateVarBasic(scip, &first.edges[i].var, name_of_variable_first_cycle.c_str(), 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY));
    SCIP_CALL(SCIPaddVar(scip, first.edges[i].var));
  
    std::string name_of_variable_second_cycle = (std::string("Sx") + std::to_string(second.edges[i].start.number) + std::to_string(first.edges[i].end.number));
    SCIP_CALL(SCIPcreateVarBasic(scip, &second.edges[i].var, name_of_variable_second_cycle.c_str(), 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY));
    SCIP_CALL(SCIPaddVar(scip, second.edges[i].var));
  }

  // Create (1) constraint
  SCIP_CALL(add_first_constraint(scip, first, second));

  // Create (2) constraint
  SCIP_CALL(add_second_constraint(scip, first, second, flags.graph_type));

  SCIP_CALL(SCIPsolve(scip));

  restoreConsole();
}