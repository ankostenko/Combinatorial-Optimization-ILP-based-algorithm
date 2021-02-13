#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string>
#include <sstream>
#include <tuple>
#include <vector>
#include <algorithm>
#include <random>

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

#include <time.h>

#include <omp.h>

#define SOKOL_TIME_IMPL
#include "sokol_time.h"

/// Chain edge fixing
/// @param ed0 Edges of degree zero
/// @param ed1 Edges of degree one
/// @param ed3 Edges of degree three
/// @param ed4 Edges of degree four
void chain_edge_fixing(Edge edge, GraphName graph_to_fix, std::vector<std::vector<Edge>> &multigraph, std::vector<int> &broken_verticies, TypeOfGraph graph_type) {
  // Fix edge
  fix_edge_in_multigraph(multigraph, broken_verticies, edge, graph_type, graph_to_fix);

  if (TypeOfGraph::DIRECTED == graph_type) {
    if (!no_three_incident_fixed_edges(broken_verticies, multigraph)) {
      return;
    }

    if (broken_verticies.empty()) {
      return;
    }

    visit_edge_in_multigraph(multigraph, edge, graph_type);

    std::vector<Edge> &outgoing_edges = multigraph[edge.start.number - 1];
    Edge outgoing_edge = { .fixed = true };
    for (Edge e: outgoing_edges) {
      if (e.fixed == false && e.start.number == edge.start.number) {
        outgoing_edge = e;
        break;
      }
    }

    if (outgoing_edge.fixed == false) {
      chain_edge_fixing(outgoing_edge, GraphName::Z_GRAPH == graph_to_fix ? GraphName::W_GRAPH : GraphName::Z_GRAPH, 
                    multigraph, broken_verticies, graph_type);
    }

    std::vector<Edge> &ingoing_edges = multigraph[edge.end.number - 1];
    Edge ingoing_edge = { .fixed = true };
    for (Edge e: ingoing_edges) {
      if (e.fixed == false && e.end.number == edge.end.number) {
        ingoing_edge = e;
        break;
      }
    }
    if (ingoing_edge.fixed == false) {
      chain_edge_fixing(ingoing_edge, GraphName::Z_GRAPH == graph_to_fix ? GraphName::W_GRAPH : GraphName::Z_GRAPH,
                    multigraph, broken_verticies, graph_type);
    }
  } else {
    // Find number of incident edges to the start vertex
    std::vector<Edge> &start_vector = multigraph[edge.start.number - 1];
    int number_of_incident_fixed_edges = 0;
    for (Edge e : start_vector) {
      if (e.graph_name == graph_to_fix) {
        number_of_incident_fixed_edges += e.fixed;
      }
    }

    if (number_of_incident_fixed_edges == 2) {
      for (Edge e: start_vector) {
        if (e.fixed) { continue; }
        chain_edge_fixing(e, GraphName::Z_GRAPH == graph_to_fix ? GraphName::W_GRAPH : GraphName::Z_GRAPH, multigraph, broken_verticies, graph_type);
      }
    }

    std::vector<Edge> &end_vector = multigraph[edge.end.number - 1];
    number_of_incident_fixed_edges = 0;
    for (Edge e: end_vector) {
      if (e.graph_name == graph_to_fix) {
        number_of_incident_fixed_edges += e.fixed;
      }
    }

    if (number_of_incident_fixed_edges == 2) {
      for (Edge e: end_vector) {
        if (e.fixed) { continue; }
        chain_edge_fixing(e, GraphName::Z_GRAPH == graph_to_fix ? GraphName::W_GRAPH : GraphName::Z_GRAPH, multigraph, broken_verticies, graph_type);
      }
    }
  }
}

/// Local search with respect to the first neighborhood
std::tuple<bool, Edge *, Edge *> local_search_1(Edge *z_graph, Edge *w_graph, Tuple<Edge> *same_edges, int attempt_limit, TypeOfGraph graph_type) {
  // Find number of components in graphs
  int number_of_components_in_z_graph = find_number_of_cycles_in_graph(z_graph, graph_type);
  int number_of_components_in_w_graph = find_number_of_cycles_in_graph(w_graph, graph_type);

  // Uncheck graphs
  for (int i = 0; i < sb_count(z_graph); i++) {
    z_graph[i].visited = false;
    w_graph[i].visited = false;
  }

  // Construct multigraph
  std::vector<std::vector<Edge>> multigraph;
  for (int i = 0; i < sb_count(z_graph); i++) {
    std::vector<Edge> vec;
    multigraph.push_back(vec);
  }
  for (int i = 0; i < sb_count(z_graph); i++) {
    std::vector<Edge> &vec1 = multigraph[(z_graph)[i].start.number - 1];
    (z_graph)[i].graph_name = GraphName::Z_GRAPH;
    vec1.push_back((z_graph)[i]);
    std::vector<Edge> &vec2 = multigraph[(z_graph)[i].end.number - 1];
    (z_graph)[i].graph_name = GraphName::Z_GRAPH;
    vec2.push_back((z_graph)[i]);
    std::vector<Edge> &vec3 = multigraph[(w_graph)[i].start.number - 1];
    (w_graph)[i].graph_name = GraphName::W_GRAPH;
    vec3.push_back((w_graph)[i]);
    std::vector<Edge> &vec4 = multigraph[(w_graph)[i].end.number - 1];
    (w_graph)[i].graph_name = GraphName::W_GRAPH;
    vec4.push_back((w_graph)[i]);
  }

  // Fix doubled edges
  fix_doubled_edges(multigraph, same_edges, graph_type);

  std::vector<int> z_graph_indices(sb_count(z_graph), 0); z_graph_indices.clear();
  for (int i = 0; i < sb_count(z_graph); i++) {
    z_graph_indices.push_back(i);
  }

  // Initialize random devices
  std::random_device rd;
  std::mt19937 g(rd());

  // Array of broken verticies
  std::vector<int> broken_verticies;

  std::shuffle(z_graph_indices.begin(), z_graph_indices.end(), g);

  Edge chosen_edge = { .start = { -1 } };
  for (int index : z_graph_indices) {
    std::vector<Edge> &z_graph_at_index = multigraph[index];
    for (Edge &e: z_graph_at_index) {
      if (e.graph_name == GraphName::Z_GRAPH && e.fixed == false && e.visited == false) {
        chosen_edge = e;
        // Mark edge as visited
        visit_edge_in_multigraph(multigraph, e, graph_type);

        // Make copy of current graph
        std::vector<std::vector<Edge>> multigraph_copy = multigraph;
        
        chain_edge_fixing(chosen_edge, GraphName::W_GRAPH, multigraph, broken_verticies, graph_type);
        
        if (graph_type == TypeOfGraph::DIRECTED && broken_verticies.empty()) {
          int new_number_of_components_in_z_graph = find_number_of_cycles_in_graph_from_multigraph(multigraph, GraphName::Z_GRAPH, graph_type);
          int new_number_of_components_in_w_graph = find_number_of_cycles_in_graph_from_multigraph(multigraph, GraphName::W_GRAPH, graph_type);

          if (new_number_of_components_in_z_graph + new_number_of_components_in_w_graph 
              < number_of_components_in_z_graph + number_of_components_in_w_graph) {
            std::tie(z_graph, w_graph) = convert_multigraph_to_two_graph(multigraph, graph_type);
            return std::make_tuple(true, z_graph, w_graph);
          } else {
            multigraph = multigraph_copy;

            // Unfix all edges and fix only doubled edges
            unfix_edges_in_multigraph(multigraph);
            fix_doubled_edges(multigraph, same_edges, graph_type);
            continue;
          }
        }
        
        std::vector<int> broken_verticies_copy = broken_verticies;
        for (int i = 0; i < attempt_limit; i++) {
          // Z contains a vertex with degree not equal to 2
          while (true) {
            // We fixed all verticies can proceed to check what solution we've got 
            if (broken_verticies.empty()) {
              break;
            }

            size_t vertex_index;
            for (size_t i = 0; i < broken_verticies.size(); i++) {
              vertex_index = broken_verticies[i];
              for (Edge e: multigraph[vertex_index - 1]) {
                if (e.fixed == false) {
                  chosen_edge = e;
                  goto found_broken_edge;
                }
              }
            }
            found_broken_edge:

            // Check if there's no three incident edges
            if (!no_three_incident_fixed_edges(broken_verticies, multigraph)) {
                goto stop_attempt;
            }

            // Check degree of vertex
            if (degree_of_vertex_in_multigraph(chosen_edge.start.number, multigraph) == 1 || 
                degree_of_vertex_in_multigraph(chosen_edge.end.number, multigraph) == 1) {
              // Degree of vertex is equal to 1

              int choose_case = rand() % 2;
              switch(choose_case) {
                case 0: {
                  chain_edge_fixing(chosen_edge, GraphName::Z_GRAPH, multigraph, broken_verticies, graph_type);
                } break;
                case 1: {
                  for (Edge e: multigraph[vertex_index - 1]) {
                    if (e.fixed == false && !are_edges_the_same_undirected_graphs(e, chosen_edge)) {
                      chain_edge_fixing(e, GraphName::Z_GRAPH, multigraph, broken_verticies, graph_type);
                      break;
                    }
                  }
                } break;
              }
            } else if (degree_of_vertex_in_multigraph(chosen_edge.start.number, multigraph) == 3 ||
                      degree_of_vertex_in_multigraph(chosen_edge.end.number, multigraph) == 3) {
              // Degree of vertex is equal to 3
              int choose_case = rand() % 2;
              switch(choose_case) {
                case 1: {
                  chain_edge_fixing(chosen_edge, GraphName::W_GRAPH, multigraph, broken_verticies, graph_type);
                } break;
                case 2: {
                  for (Edge e: multigraph[vertex_index - 1]) {
                    if (e.fixed == false && !are_edges_the_same_undirected_graphs(e, chosen_edge)) {
                      chain_edge_fixing(e, GraphName::W_GRAPH, multigraph, broken_verticies, graph_type);
                      break;
                    }
                  }
                }
              }
            }

            if (degree_of_vertex_in_multigraph(chosen_edge.start.number, multigraph) == 4 ||
                degree_of_vertex_in_multigraph(chosen_edge.end.number, multigraph) == 4) {
              return std::make_tuple(false, z_graph, w_graph);
            } 
            
            if (degree_of_vertex_in_multigraph(chosen_edge.start.number, multigraph) == 0 || 
                degree_of_vertex_in_multigraph(chosen_edge.end.number, multigraph) == 0) {
              return std::make_tuple(false, z_graph, w_graph);
            }

            set_brokeness_of_verticies(multigraph[chosen_edge.start.number - 1], chosen_edge.start.number, broken_verticies);
            set_brokeness_of_verticies(multigraph[chosen_edge.end.number - 1], chosen_edge.end.number, broken_verticies);
          }

          // If number of components in the graphs has decreased stop attempts
          {
            int new_number_of_components_in_z_graph = 0;
            int new_number_of_components_in_w_graph = 0;
            {
                new_number_of_components_in_z_graph = find_number_of_cycles_in_graph_from_multigraph(multigraph, GraphName::Z_GRAPH, graph_type);
                new_number_of_components_in_w_graph = find_number_of_cycles_in_graph_from_multigraph(multigraph, GraphName::W_GRAPH, graph_type);
            }
            if (new_number_of_components_in_z_graph + new_number_of_components_in_w_graph 
                < number_of_components_in_z_graph + number_of_components_in_w_graph) {
              std::tie(z_graph, w_graph) = convert_multigraph_to_two_graph(multigraph, graph_type);
              return std::make_tuple(true, z_graph, w_graph);
            }
          }

          stop_attempt:
          // Restore original Z and W graphs
          multigraph = multigraph_copy;
          broken_verticies = broken_verticies_copy;

          // Unfix all edges and fix only doubled edges
          unfix_edges_in_multigraph(multigraph);
          fix_doubled_edges(multigraph, same_edges, graph_type);
        }
      }
    }
  }

  return std::make_tuple(false, z_graph, w_graph);
}

/// Variable neighborhood descent
std::tuple<bool, Edge*, Edge*> variable_neighborhood_descent(SCIP *scip, Edge *z_graph, Edge *w_graph, Tuple<Edge> *same_edges, int attempt_limit, TypeOfGraph graph_type) {
  do {
    // Local search
    bool has_improvement = false;
    std::tie(has_improvement, z_graph, w_graph) = local_search_1(z_graph, w_graph, same_edges, attempt_limit, graph_type);
    
    // Determine if Z graph and W graph are hamilton decomposition
    if (cycles_are_new_decomposition(scip, z_graph, w_graph, same_edges, graph_type)) {
      return std::make_tuple(true, z_graph, w_graph);
    }
    
    if (has_improvement) {
      // Improvement with respect to the first neighborhood
      continue;
    }
    break;
  } while(true);

  return std::make_tuple(false, z_graph, w_graph);
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
  double local_search_1_analysis_average_time = 0.0;
  double has_solution_average_time = 0.0;
  double no_solution_average_time = 0.0;
  float number_of_iterations = 0.0f;

  int number_of_tests = 0;

  std::vector<int> vertex_vector;
  vertex_vector.reserve(flags.number_of_verticies); vertex_vector.clear();
  for (int i = 1; i <= flags.number_of_verticies; i++) {
    vertex_vector.push_back(i);
  }

  // Initialize random devices
  std::random_device rd;
  std::mt19937 g(rd());

  while(true) {
    BundleOfGraphs bundle;

    if (flags.generate_cycles) {
      std::shuffle(vertex_vector.begin(), vertex_vector.end(), g);
      for (int vertex: vertex_vector) {
        bundle.first.add_vertex(Vertex({vertex}));
      }

      bundle.first = create_edges_from_verticies(&bundle.first);

      std::shuffle(vertex_vector.begin(), vertex_vector.end(), g);
      for (int vertex: vertex_vector) {
        bundle.second.add_vertex(Vertex({vertex}));
      }
      bundle.second = create_edges_from_verticies(&bundle.second);
    } else {
      // Load cycles to the program
      if (read_cycles_from_stream(&bundle, file, total_file_size) == false) {
        break;
      }
    }

    uint64_t general_start =  stm_now();
    uint64_t has_solution_start = stm_now();
    uint64_t no_solution_start = stm_now();

    // Increment number of tests
    number_of_tests += 1;

    // Check if number of tests with generated cycles is got their limit
    if (flags.generate_cycles) {
      if (number_of_tests > flags.number_of_tests) {
        number_of_tests -= 1;
        break;
      }
    }

    for (int i = number_of_tests; i > 0; i /= 2) {
      printf("\b");
    }
    printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
    printf("Running test #%d", number_of_tests);

    // Set graphs to the corresponding cycles
    Graph first = bundle.first;
    Graph second = bundle.second;

    // Start constraint construction timings
    uint64_t constraint_construction_time_start = stm_now();

    SCIP *scip = NULL;

    // Initialize SCIP
    init_scip(&scip);

    // Set messages to quite
    SCIPsetMessagehdlrQuiet(scip, true);

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

    // Stop constraint construction timing
    constraint_construction_time_average_time += stop_time(constraint_construction_time_start);

    do {
      // Number of iterations of the algorithm for one cycle 
      number_of_iterations += 1.0f;

      Edge *z_graph = NULL;
      Edge *w_graph = NULL;

      // Start scip timing
      uint64_t scip_start = stm_now();

      // Solve problem
      SCIP_CALL(SCIPsolve(scip));

      // Stop scip timing
      scip_average_time += stop_time(scip_start);
    
      // Print solution on the current step
      if (!(SCIPgetNSols(scip) > 0)) {
        // Stop general timer and no solution timer
        no_solution_average_time += stop_time(no_solution_start);

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

      // Start cycle analysis time
      uint64_t cycle_analysis_start = stm_now();

      if (cycles_are_new_decomposition(scip, z_graph, w_graph, same_edges, flags.graph_type)) {
        // Found new decomposition
        has_decomposition++;

        has_solution_average_time += stop_time(has_solution_start);
        cycle_analysis_average_time += stop_time(cycle_analysis_start);
        break;
      }
      
      // Stop cycle analysis timer
      cycle_analysis_average_time += stop_time(cycle_analysis_start);

      // Start local search time
      uint64_t local_search_1_analysis_start = stm_now();

      // Local search with respect to the first neighborhood
      if (flags.first_neighborhood_enabled) {
        bool has_solution = false;
        std::tie(has_solution, z_graph, w_graph) = variable_neighborhood_descent(scip, z_graph, w_graph, same_edges, flags.attempt_limit, flags.graph_type);
        
        bool z_graph_is_different = false;
        bool w_graph_is_different = false;
        if (has_solution) {
          for (int i = 0; i < sb_count(z_graph); i++) {
            if (!contains_edge(bundle.first.edges[i], z_graph, flags.graph_type)) {
              z_graph_is_different = true;
              break;
            }
          }

          for (int i = 0; i < sb_count(w_graph); i++) {
            if (!contains_edge(bundle.first.edges[i], w_graph, flags.graph_type)) {
              w_graph_is_different = true;
              break;
            }
          }
        }

        if (z_graph_is_different && w_graph_is_different) {
          // We have a valid solution
          local_search_1_analysis_average_time += stop_time(local_search_1_analysis_start);
          general_average_time += stop_time(general_start);
          has_solution_average_time += stop_time(has_solution_start);
          
          has_decomposition++;
          goto has_valid_solution;
        }
      }

      // Stop local search timer
      local_search_1_analysis_average_time += stop_time(local_search_1_analysis_start);

      // Free z and w graphs
      sb_free(z_graph);
      sb_free(w_graph);

    } while(true);

    general_average_time += stop_time(general_start);
    has_valid_solution:

    // Release variables
    for (size_t i = 0; i < bundle.first.number_of_edges(); i++) {
      SCIP_CALL(SCIPreleaseVar(scip, &bundle.first.edges[i].var));
      SCIP_CALL(SCIPreleaseVar(scip, &bundle.second.edges[i].var));
    }

    // Free cycles
    sb_free(bundle.first.edges);
    sb_free(bundle.first.verticies);
    sb_free(bundle.second.edges);
    sb_free(bundle.second.verticies);

    // Release SCIP structure
    SCIP_CALL(SCIPfree(&scip));
  }

  printf("\nResults:\nNo decomposition: %d, Decomposition: %d\n", no_decomposition, has_decomposition);
  printf("Number of tests: %d\n", number_of_tests);
  printf("General average time: %f ms\n", general_average_time / number_of_tests);
  printf("Has solution average time: %f ms\n", has_solution_average_time / has_decomposition);
  printf("No solution average time: %f ms\n", no_solution_average_time / no_decomposition);
  printf("Constraint construction average time: %f ms\n", constraint_construction_time_average_time / number_of_tests);
  printf("SCIP average time: %f ms\n", scip_average_time / number_of_tests);
  printf("Cycle analysis average time: %f ms\n", cycle_analysis_average_time / number_of_tests);
  printf("Local search 1 average time: %f ms\n", local_search_1_analysis_average_time / number_of_tests);
  printf("Average number of iterations: %f\n", number_of_iterations / number_of_tests);

  restoreConsole();
}