#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string>
#include <sstream>
#include <tuple>

#include "utils.h"
#include "config_param_parsing.h"
#include "win_console_color.h"
#include "graph.h"
#include "cycle_extraction.h"
#include "print/printing.h"

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/struct_prob.h"
#include "scip/struct_scip.h"
#include "scip/sol.h"

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
  SCIPinfoMessage(scip, NULL,  "\n");

  SCIP_CALL(SCIPreleaseCons(scip, &constraint));

  return SCIP_OKAY;
}

/// Add second constraint to the problem
SCIP_RETCODE add_second_constraint(SCIP *scip, Graph first, Graph second, TypeOfGraph graph_type) {
  if (graph_type == TypeOfGraph::DIRECTED) {
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
          SCIPinfoMessage(scip, NULL,  "\n");


          SCIP_CALL(SCIPreleaseCons(scip, &constraint_ingoing_edges));

          // Outgoing edges
          SCIP_CALL(SCIPcreateConsBasicLinear(scip, &constraint_outgoing_edges, "Ingoing edges", 2, variables_outgoing_edges, coefficients, 1.0, 1.0));
          SCIP_CALL(SCIPaddCons(scip, constraint_outgoing_edges));
          
          // TODO: remove this print
          SCIP_CALL(SCIPprintCons(scip, constraint_outgoing_edges, NULL));
          SCIPinfoMessage(scip, NULL,  "\n");

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
          SCIPinfoMessage(scip, NULL,  "\n");

          SCIP_CALL(SCIPreleaseCons(scip, &constraint));

          // Move to the next vertex
          break;
        }
      }
    }
  }
  
  return SCIP_OKAY;
}

/// Check whether edge contains inside the edges array
/// @param edges Array of edges to look in
/// @param edge Edge to compare with
/// @param graph_type Type of graph
/// @return Returns true if edge contains in the array, otherwise returns false
bool is_edge_contains_in(Tuple<Edge> *edges, Edge edge, TypeOfGraph graph_type) {
  if (graph_type == TypeOfGraph::DIRECTED) {
    for (int i = 0; i < sb_count(edges); i++) {
      if (are_edges_the_same_directed_graphs(edges[i].first, edge)) {
        return true;
      }
    }
  } else if (graph_type == TypeOfGraph::UNDIRECTED) {
    for (int i = 0; i < sb_count(edges); i++) {
      if (are_edges_the_same_undirected_graphs(edges[i].first, edge)) {
        return true;
      }
    }
  } else {
    assert(0 && !"Unknown type of graph");
  }

  return false;
}

/// Add (3) or (4) constraint to the problem
/// @param scip SCIP problem formulation variable
/// @param graph Cycle which is used to construct constraint
/// @param same_edges Array of edges that are present in both cycles
/// @param graph_type Type of graph
/// @return SCIP_RETCODE
SCIP_RETCODE add_third_or_fourth_constraint(SCIP *scip, Graph graph, Tuple<Edge> *same_edges, TypeOfGraph graph_type, const char *name_of_constraint) {
  SCIP_CONS *constraint = NULL;
  SCIP_VAR **constraint_variables = NULL;

  int number_of_verticies = graph.number_of_verticies();

  for (size_t i = 0; i < graph.number_of_edges(); i++) {
    Edge edge = graph.edges[i];
    // If edge is not in the same_edges array
    if (is_edge_contains_in(same_edges, edge, graph_type) == false) {
      sb_push(constraint_variables, edge.var);
    }
  }

  assert((number_of_verticies - sb_count(same_edges)) == sb_count(constraint_variables));

  double *values = { };
  for (int i = 0; i < sb_count(constraint_variables); i++) {
    sb_push(values, 1.0);
  }

  SCIP_CALL(SCIPcreateConsBasicLinear(scip, &constraint, name_of_constraint, 
      sb_count(constraint_variables), constraint_variables, values, -SCIPinfinity(scip), number_of_verticies - sb_count(same_edges) - 2));
  SCIP_CALL(SCIPaddCons(scip, constraint));
  
  // TODO: remove
  SCIP_CALL(SCIPprintCons(scip, constraint, NULL));
  SCIPinfoMessage(scip, NULL,  "\n");

  SCIP_CALL(SCIPreleaseCons(scip, &constraint));

  return SCIP_OKAY;
}

/// Find same edges
/// @param first First cycle
/// @param second Second cycle
/// @param graph_type Type of graph
/// @return Array of same edges for both cycles
Tuple<Edge> *find_same_edges(Graph first, Graph second, TypeOfGraph graph_type) {
  Tuple<Edge> *same_edges = NULL;
  
  if (graph_type == TypeOfGraph::DIRECTED) {
    for (size_t i = 0; i < first.number_of_edges(); i++) {
      Edge first_cycle_edge = first.edges[i];
      for (size_t j = 0; j < second.number_of_edges(); j++) {
        Edge second_cycle_edge = second.edges[j];

        // Check edges for equality
        if (are_edges_the_same_directed_graphs(first_cycle_edge, second_cycle_edge)) {
          sb_push(same_edges, Tuple<Edge>({first_cycle_edge, second_cycle_edge}));
          break;
        }
      }
    }
  } else if (graph_type == TypeOfGraph::UNDIRECTED) {
    for (size_t i = 0; i < first.number_of_edges(); i++) {
      Edge first_cycle_edge = first.edges[i];
      for (size_t j = 0; j < second.number_of_edges(); j++) {
        Edge second_cycle_edge = second.edges[j];

        // Check edges for equality
        if (are_edges_the_same_undirected_graphs(first_cycle_edge, second_cycle_edge)) {
          sb_push(same_edges, Tuple<Edge>({first_cycle_edge, second_cycle_edge}));
          break;
        }
      }
    }
  } else {
    assert(0 && !"Unknown type of graph");
  }

  return same_edges;
}

/// (7) Constraint
/// Separate same edges to different cycles
SCIP_RETCODE add_seventh_constraint(SCIP *scip, Tuple<Edge> *same_edges) {
  SCIP_CONS *constraint;
  SCIP_Real coefficient = 1.0;

  for (int i = 0; i < sb_count(same_edges); i++) {
    Tuple<Edge> tp = same_edges[i];

    Edge z_cycle_edge = tp.first;
    SCIP_CALL(SCIPcreateConsBasicLinear(scip, &constraint, "Separation Constraint. Z cycle", 1, &z_cycle_edge.var, &coefficient, 1.0, 1.0));
    SCIP_CALL(SCIPaddCons(scip, constraint));

    // TODO: remove
    SCIP_CALL(SCIPprintCons(scip, constraint, NULL));
    SCIPinfoMessage(scip, NULL,  "\n");

    SCIP_CALL(SCIPreleaseCons(scip, &constraint));

    Edge w_cycle_edge = tp.second;
    SCIP_CALL(SCIPcreateConsBasicLinear(scip, &constraint, "Separation Constraint. W cycle", 1, &w_cycle_edge.var, &coefficient, 0.0, 0.0));
    SCIP_CALL(SCIPaddCons(scip, constraint));
    
    // TODO: remove
    SCIP_CALL(SCIPprintCons(scip, constraint, NULL));
    SCIPinfoMessage(scip, NULL,  "\n");

    SCIP_CALL(SCIPreleaseCons(scip, &constraint));
  }

  return SCIP_OKAY;
}

/// Find cycle
/// Returns newly allocated cycle.
Edge* find_cycle(Edge *graph, TypeOfGraph graph_type) {
  Edge *cycle = NULL;
  
  if (graph_type == TypeOfGraph::UNDIRECTED) {
    // Branch for undirected graphs
    for (int i = 0; i < sb_count(graph); i++) {
      Edge current_edge = graph[i];
      if (current_edge.visited == true) { continue; }
      current_edge.visited = true;

      for (int j = i; j < sb_count(graph); j++) {
        int edge_index = (j + 1 == sb_count(graph)) ? 0 : j + 1;
        Edge next_edge = graph[edge_index];
  
        if (next_edge.visited == true) { continue; }

        if (edges_conjuncted_directed_graphs(current_edge, next_edge)) {
          // Edges are conjuncted so we add them to the cycle
          sb_push(cycle, next_edge);
          current_edge = next_edge;
          current_edge.visited = true;
          j = i;
        }    
      }
    }
  } else if (graph_type == TypeOfGraph::DIRECTED) {
    // Branch for directed graphs
        // Branch for undirected graphs
    for (int i = 0; i < sb_count(graph); i++) {
      Edge current_edge = graph[i];
      if (current_edge.visited == true) { continue; }
      current_edge.visited = true;

      for (int j = i; j < sb_count(graph); j++) {
        int edge_index = (j + 1 == sb_count(graph)) ? 0 : j + 1;
        Edge next_edge = graph[edge_index];
  
        if (next_edge.visited == true) { continue; }

        if (edges_conjuncted_undirected_graphs(current_edge, next_edge)) {
          // Edges are conjuncted so we add them to the cycle
          sb_push(cycle, next_edge);
          current_edge = next_edge;
          current_edge.visited = true;
          j = i;
        }    
      }
    }
  } else {
    assert(0 && !"We have only two types of graphs");
  }

  return cycle;
}

/// Find components of a graph and add them to constraints of the problem
/// @param graph Graph to find 
/// @return Return true if found cycle is hamilton cycle of the graph, false otherwise
bool find_cycles_and_add_to_constraints(SCIP *scip, Edge *graph, Tuple<Edge> *same_edges, TypeOfGraph graph_type) {
  // Try to find cycles until all edges are not available
  while (true) {
    Edge *current_cycle = find_cycle(graph, graph_type);

    if (current_cycle == NULL) {
      return false;
    }

    // Iterate over the current cycle and add this cycle
    // to the constraints
    SCIP_CONS *constraint = NULL;
    SCIP_VAR *variables = NULL;    

    // Current cycle is Hamilton cycle
    if (sb_count(current_cycle) == sb_count(graph)) {
      return true;
    }

    // Free current cycle
    sb_free(current_cycle);
  }

  return false;
}

/// Analyzes Z-graph and W-graph 
/// @return If Z-graph and W-graph are new Hamilton Decomposition return true otherwise false
bool cycles_are_new_decomposition(SCIP *scip, Edge *z_graph, Edge *w_graph, Tuple<Edge> *same_edges, TypeOfGraph graph_type) {
  bool z_graph_is_complete_cycle = find_cycles_and_add_to_constraints(scip, z_graph, same_edges, graph_type);
  bool w_graph_is_complete_cycle = find_cycles_and_add_to_constraints(scip, w_graph, same_edges, graph_type);
  
  return z_graph_is_complete_cycle && w_graph_is_complete_cycle;
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

  const char *path_to_the_file = "test_cases\\test_6_undir_paper_motor.txt";

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

  Edge *z_graph = NULL;
  Edge *w_graph = NULL;

  do {
    // Solve problem
    SCIP_CALL(SCIPsolve(scip));
  
    // Print solution on the current step
    if (SCIPgetNSols(scip) > 0) {
      SCIPinfoMessage(scip, NULL, "\nSolution:\n");
      SCIP_CALL(SCIPprintSol(scip, SCIPgetBestSol(scip), NULL, TRUE));
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

    if (cycles_are_new_decomposition(scip, z_graph, w_graph, same_edges, flags.graph_type)) {
      // Found new decomposition
      break;
    }

    sb_free(z_graph);
    sb_free(w_graph);
  } while(true);

  restoreConsole();
}