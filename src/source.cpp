#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string>
#include <sstream>
#include <tuple>

#include "stretchy_buffer.h"

#include "utils.h"
#include "config_param_parsing.h"
#include "win_console_color.h"
#include "graph.h"
#include "cycle_extraction.h"
#include "print/printing.h"
#include "ilp.h"

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/struct_prob.h"
#include "scip/struct_scip.h"
#include "scip/sol.h"

#define SOKOL_TIME_IMPL
#include "sokol_time.h"

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

  const char *path_to_the_file = flags.path_to_test_file;

  FILE *file = NULL;
  int error_code = fopen_s(&file, path_to_the_file, "r");
  if (error_code != 0) {
    printf("Path to the file: %s\n", path_to_the_file);
    report_problem_to_user(true, "Couldn't open the file with cycles");
  }
  size_t total_file_size = get_size_of_file(file);

  int no_decomposition = 0;
  int has_decomposition = 0;

  // Init sokol time
  stm_setup();

  double general_average_time = 0.0;
  double constraint_construction_time_average_time = 0.0;
  double scip_average_time = 0.0;
  double cycle_analysis_average_time = 0.0;
  float number_of_iterations = 0.0f;

  int number_of_tests = 0;

  while(true) {
    uint64_t general_start =  stm_now();

    // Load cycles to the program
    BundleOfGraphs bundle;
    if (read_cycles_from_stream(&bundle, file, total_file_size) == false) {
      break;
    }

    // Increment number of tests
    number_of_tests += 1;

    // Set graphs to the corresponding cycles
    Graph first = bundle.first;
    Graph second = bundle.second;

    uint64_t constraint_construction_time_start = stm_now();

    SCIP *scip = NULL;

    // Initialize SCIP
    init_scip(&scip);

    // Add variables to the problem from cycles
    Edge *all_edges = NULL;
    for (size_t i = 0; i < first.number_of_verticies(); i++) {
      std::string name_of_variable_first_cycle = (std::string("Fx") + std::to_string(first.edges[i].start.number) + std::to_string(first.edges[i].end.number));
      SCIP_CALL(SCIPcreateVarBasic(scip, &first.edges[i].var, name_of_variable_first_cycle.c_str(), 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY));
      SCIP_CALL(SCIPaddVar(scip, first.edges[i].var));

      // Add first edge to the array of all edges
      sb_push(all_edges, first.edges[i]);

      std::string name_of_variable_second_cycle = (std::string("Sx") + std::to_string(second.edges[i].start.number) + std::to_string(second.edges[i].end.number));
      SCIP_CALL(SCIPcreateVarBasic(scip, &second.edges[i].var, name_of_variable_second_cycle.c_str(), 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY));
      SCIP_CALL(SCIPaddVar(scip, second.edges[i].var));
    
      // Add second edge to the array of all edges
      sb_push(all_edges, second.edges[i]);
    }

    // Add (1) constraint
    SCIP_CALL(add_first_constraint(scip, first, second));

    // Add (2) constraint
    SCIP_CALL(add_second_constraint(scip, first, second, flags.graph_type));

    // Find the same edges of the cycles
    Tuple<Edge> *same_edges = find_same_edges(first, second, flags.graph_type);

    // Add (3) constraint
    SCIP_CALL(add_third_or_fourth_constraint(scip, first, same_edges, flags.graph_type, "Third constraint"));

    // Add (4) constraint
    SCIP_CALL(add_third_or_fourth_constraint(scip, second, same_edges, flags.graph_type, "Fourth constraint"));
    
    // Add (7) constraint
    SCIP_CALL(add_seventh_constraint(scip, same_edges));

    uint64_t constraint_construction_time_elapsed = stm_since(constraint_construction_time_start);
    double constraint_construction_milliseconds = stm_ms(constraint_construction_time_elapsed);

    constraint_construction_time_average_time += constraint_construction_milliseconds;

    do {
      // Number of iterations of the algorithm for one cycle 
      number_of_iterations += 1.0f;

      Edge *z_graph = NULL;
      Edge *w_graph = NULL;

      uint64_t scip_start = stm_now();

      // Solve problem
      SCIP_CALL(SCIPsolve(scip));

      uint64_t scip_elapsed = stm_since(scip_start);
      double scip_milliseconds = stm_ms(scip_elapsed);
      scip_average_time += scip_milliseconds;
    
      // Print solution on the current step
      if (!(SCIPgetNSols(scip) > 0)) {
        // No solution
        no_decomposition++;
        break;
      }

      // Analyze solution vector and separate edges to Z and W cycles
      SCIP_PROB *original_problem = scip->origprob;
      for (int i = 0; i < original_problem->nvars; i++) {
        SCIP_Real variable_value = SCIPsolGetVal(SCIPgetBestSol(scip), scip->set, scip->stat, original_problem->vars[i]);

        Edge edge = all_edges[i];
        if (fabs(variable_value - 1.0) < EPSILON) {
          // Variable value is 1.0 so the variable goes to Z cycle
          sb_push(z_graph, edge);
        } else {
          // Variable value is 0.0 so the variable goes to W cycle
          sb_push(w_graph, edge);
        }
      }

      // Free transformed problem in order to add new constraint to the problem
      SCIPfreeTransform(scip);

      uint64_t cycle_analysis_start = stm_now();

      if (cycles_are_new_decomposition(scip, z_graph, w_graph, same_edges, flags.graph_type)) {
        // Found new decomposition
        has_decomposition++;
        
        uint64_t cycle_analysis_elapsed = stm_since(cycle_analysis_start);
        double cycle_analysis_milliseconds = stm_ms(cycle_analysis_elapsed);
        cycle_analysis_average_time += cycle_analysis_milliseconds;

        break;
      }

      uint64_t cycle_analysis_elapsed = stm_since(cycle_analysis_start);
      double cycle_analysis_milliseconds = stm_ms(cycle_analysis_elapsed);
      cycle_analysis_average_time += cycle_analysis_milliseconds;

      sb_free(z_graph);
      sb_free(w_graph);

      uint64_t general_elapsed = stm_since(general_start);
      double general_milliseconds = stm_ms(general_elapsed);

      general_average_time += general_milliseconds;

    } while(true);
  }

  printf("Results:\nNo decomposition: %d, Decomposition: %d\n", no_decomposition, has_decomposition);
  printf("General average time: %f ms\n", general_average_time / number_of_tests);
  printf("Constraint construction average time: %f ms\n", constraint_construction_time_average_time / number_of_tests);
  printf("SCIP average time: %f ms\n", scip_average_time / number_of_tests);
  printf("Cycle analysis average time: %f ms\n", cycle_analysis_average_time / number_of_tests);
  printf("Average number of iterations: %f\n", number_of_iterations / number_of_tests);

  restoreConsole();
}