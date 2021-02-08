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

#define SOKOL_TIME_IMPL
#include "sokol_time.h"


void set_brokeness_of_verticies(std::vector<Edge> edges, int vertex_number, std::vector<int> &broken_verticies) {
  // Check if the verticies now broken
  // 1. Find degree of the vertex
  int degree_of_vertex = 0;
  for (Edge e: edges) {
    if (e.graph_name == GraphName::Z_GRAPH) {
      degree_of_vertex += 1;
    }
  }

  if (degree_of_vertex == 2) {
    // 2. If degree equal to 2 remove vertex from broken verticies
    auto it = std::find(broken_verticies.begin(), broken_verticies.end(), vertex_number);
    if (it != broken_verticies.end()) {
      broken_verticies.erase(it);
    }
  } else {
    // 3. If degree not equal to 2 add vertex if it not there already
    auto it = std::find(broken_verticies.begin(), broken_verticies.end(), vertex_number);
    if (it == broken_verticies.end()) {
      broken_verticies.push_back(vertex_number);
    }
  }
}

/// Fix edge in multigraph and set in which graph it goes
void fix_edge_in_multigraph(std::vector<std::vector<Edge>> &multigraph, std::vector<int> &broken_verticies, Edge edge, 
                                                                            TypeOfGraph graph_type, GraphName graph_name) {
  std::vector<Edge> start_vector = multigraph[edge.start.number - 1];
  std::vector<Edge> end_vector = multigraph[edge.end.number - 1];

  for (Edge &e : start_vector) {
    if (TypeOfGraph::DIRECTED == graph_type) {
      if (are_edges_the_same_directed_graphs(e, edge)) {
        e.fixed = true;
        e.graph_name = graph_name;
      }
    } else {
      if (are_edges_the_same_undirected_graphs(e, edge)) {
        e.fixed = true;
        e.graph_name = graph_name;
      }
    }
  }

  for (Edge &e: end_vector) {
    if (TypeOfGraph::DIRECTED == graph_type) {
      if (are_edges_the_same_directed_graphs(e, edge)) {
        e.fixed = true;
        e.graph_name = graph_name;
      }
    } else {
      if (are_edges_the_same_undirected_graphs(e, edge)) {
        e.fixed = true;
        e.graph_name = graph_name;
      }
    }
  }

  set_brokeness_of_verticies(start_vector, edge.start.number, broken_verticies);
  set_brokeness_of_verticies(end_vector, edge.end.number, broken_verticies);
}

/// Chain edge fixing
/// @param ed0 Edges of degree zero
/// @param ed1 Edges of degree one
/// @param ed3 Edges of degree three
/// @param ed4 Edges of degree four
void chain_edge_fixing(Edge edge, GraphName graph_to_fix, std::vector<std::vector<Edge>> &multigraph, std::vector<int> &broken_verticies, TypeOfGraph graph_type) {
  // Fix edge
  fix_edge_in_multigraph(multigraph, broken_verticies, edge, graph_type, graph_to_fix);

  if (TypeOfGraph::DIRECTED) {

  } else {
    // Find number of incident edges to the start vertex
    std::vector<Edge> start_vector = multigraph[edge.start.number - 1];
    int number_of_incident_fixed_edges = 0;
    for (Edge e : start_vector) {
      number_of_incident_fixed_edges += e.fixed;
    }

    if (number_of_incident_fixed_edges == 2) {
      for (Edge e: start_vector) {
        if (e.fixed) { continue; }
        chain_edge_fixing(e, GraphName::Z_GRAPH == graph_to_fix ? GraphName::W_GRAPH : GraphName::Z_GRAPH, multigraph, broken_verticies, graph_type);
      }
    }

    std::vector<Edge> end_vector = multigraph[edge.end.number - 1];
    number_of_incident_fixed_edges = 0;
    for (Edge e: end_vector) {
      number_of_incident_fixed_edges += e.fixed;
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
bool local_search_1(Edge *z_graph, Edge *w_graph, Tuple<Edge> *same_edges, int attempt_limit, TypeOfGraph graph_type) {
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
    std::vector<Edge> &vec1 = multigraph[z_graph[i].start.number - 1];
    z_graph[i].graph_name = GraphName::Z_GRAPH;
    vec1.push_back(z_graph[i]);
    std::vector<Edge> &vec2 = multigraph[z_graph[i].end.number - 1];
    z_graph[i].graph_name = GraphName::Z_GRAPH;
    vec2.push_back(z_graph[i]);
    std::vector<Edge> &vec3 = multigraph[w_graph[i].start.number - 1];
    w_graph[i].graph_name = GraphName::W_GRAPH;
    vec3.push_back(w_graph[i]);
    std::vector<Edge> &vec4 = multigraph[w_graph[i].end.number - 1];
    w_graph[i].graph_name = GraphName::W_GRAPH;
    vec4.push_back(w_graph[i]);
  }

  // Fix doubled edges
  for (auto &vec: multigraph) {
    for (Edge &e: vec) {
      if (is_edge_contains_in(same_edges, e, graph_type)) {
        e.fixed = true;
      }
    }
  }

  for (auto vec : multigraph) {
    for (Edge edge : vec) {
      printf("%d-%d ", edge.start.number, edge.end.number);
    }
    printf("\n");
  }

  std::vector<int> z_graph_indices(sb_count(z_graph), 0); z_graph_indices.clear();
  for (int i = 0; i < sb_count(z_graph); i++) {
    z_graph_indices.push_back(i);
  }

  std::random_device rd;
  std::mt19937 g(rd());

  std::vector<int> broken_verticies;

  while (true) {
    // Find unchecked and unfixed edge of z
    // 1. Shuffle indices
    std::shuffle(z_graph_indices.begin(), z_graph_indices.end(), g);
    Edge chosen_edge = { .start = { -1 } };
    // 2. Choose unfixed and unchecked edge
    for (int index : z_graph_indices) {
      std::vector<Edge> z_graph_at_index = multigraph[index];
      for (Edge e: z_graph_at_index) {
        if (e.fixed == false && e.visited == false) {
          chosen_edge = e;
        }
      }
    }
    // 3. Break if we don't have unchecked and unfixed edges
    if (chosen_edge.start.number == -1) { break; }

    chain_edge_fixing(chosen_edge, GraphName::Z_GRAPH, multigraph, broken_verticies, graph_type);
    for (int i = 0; i < attempt_limit; i++) {
      while (true) {
        // We fixed all verticies can proceed to check what solution we got 
        if (broken_verticies.empty()) {
          break;
        }

        // Check degree of vertex
        if (degree_of_vertex_in_multigraph(chosen_edge.start.number, multigraph) == 1) {
          // Degree of vertex is equal to 1
          switch(1) {
            case 1: {
              std::vector<Edge> edges = multigraph[chosen_edge.start.number - 1];
              std::shuffle(edges.begin(), edges.end(), g);
              for (Edge e: edges) {
                if (e.fixed == false) {
                  chain_edge_fixing(e, GraphName::W_GRAPH, multigraph, broken_verticies, graph_type);
                }
              }
            } break;
          }
        }
      }
    }
  }

  return false;
}

// bool local_search_2(Edge *z_graph, Edge *w_graph, int depth_limit, TypeOfGraph graph_type) {
//   return false;
// }

/// Variable neighborhood descent
void variable_neighborhood_descent(Edge *z_graph, Edge *w_graph, Tuple<Edge> *same_edges, int attempt_limit, int depth_limit, TypeOfGraph graph_type) {
  do {
    // Determine if Z graph and W graph are hamilton decomposition
    Edge *z_cycle = find_cycle(z_graph, graph_type);
    if (sb_count(z_cycle) == sb_count(z_graph)) {
      Edge *w_cycle = find_cycle(w_graph, graph_type);
      if (sb_count(w_cycle) == sb_count(w_graph)) {
        sb_free(z_cycle);
        sb_free(w_cycle);
        return;
      }
    }

    // Local search
    if (local_search_1(z_graph, w_graph, same_edges, attempt_limit, graph_type)) {
      // Improvement with respect to the first neighborhood
      continue;
    }

    // local_search_2(z_graph, w_graph, depth_limit, graph_type);
    // Z and W is a local minimum with respect to both neighborhood structures
    break;
  } while(true);
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
  Edge *z_graph = NULL;
  Edge *w_graph = NULL;
  
  Edge edge = { .start = { 1 }, .end = { 5 } };
  sb_push(z_graph, edge);
  edge = { .start = { 5 }, .end = { 8 } };
  sb_push(z_graph, edge);
  edge = { .start = { 8 }, .end = { 4 } };
  sb_push(z_graph, edge);
  edge = { .start = { 4 }, .end = { 1 } };
  sb_push(z_graph, edge);
  edge = { .start = { 6 }, .end = { 2 } };
  sb_push(z_graph, edge);
  edge = { .start = { 2 }, .end = { 3 } };
  sb_push(z_graph, edge);
  edge = { .start = { 3 }, .end = { 7 } };
  sb_push(z_graph, edge);
  edge = { .start = { 7 }, .end = { 6 } };
  sb_push(z_graph, edge);

  edge = { .start = { 1 }, .end = { 5 } };
  sb_push(w_graph, edge);
  edge = { .start = { 5 }, .end = { 6 } };
  sb_push(w_graph, edge);
  edge = { .start = { 6 }, .end = { 2 } };
  sb_push(w_graph, edge);
  edge = { .start = { 2 }, .end = { 1 } };
  sb_push(w_graph, edge);
  edge = { .start = { 4 }, .end = { 8 } };
  sb_push(w_graph, edge);
  edge = { .start = { 8 }, .end = { 7 } };
  sb_push(w_graph, edge);
  edge = { .start = { 7 }, .end = { 3 } };
  sb_push(w_graph, edge);
  edge = { .start = { 3 }, .end = { 4 } };
  sb_push(w_graph, edge);

  Graph first  = { .edges = z_graph };
  Graph second = { .edges = w_graph };

  Tuple<Edge> *same_edges = find_same_edges(first, second, TypeOfGraph::UNDIRECTED);

  local_search_1(z_graph, w_graph, same_edges, 5, TypeOfGraph::UNDIRECTED);

  return 0;


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

      printf("Z\n");
      print_edges(z_graph);
      printf("W\n");
      print_edges(w_graph);

      // Local search with respect to the first neighborhood
      local_search_1(z_graph, w_graph, same_edges, 1, flags.graph_type);

      // Free z and w graphs
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