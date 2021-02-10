#include "ilp.h"

#include "stretchy_buffer.h"

#include "print/printing.h"
#include "utils.h"

#include "scip/scip.h"
#include "scip/cons_linear.h"
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

          SCIP_CALL(SCIPreleaseCons(scip, &constraint_ingoing_edges));

          // Outgoing edges
          SCIP_CALL(SCIPcreateConsBasicLinear(scip, &constraint_outgoing_edges, "Outgoing edges", 2, variables_outgoing_edges, coefficients, 1.0, 1.0));
          SCIP_CALL(SCIPaddCons(scip, constraint_outgoing_edges));
          
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

          SCIP_CALL(SCIPreleaseCons(scip, &constraint));

          // Move to the next vertex
          break;
        }
      }
    }
  }
  
  return SCIP_OKAY;
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
  
  SCIP_CALL(SCIPreleaseCons(scip, &constraint));

  return SCIP_OKAY;
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

    SCIP_CALL(SCIPreleaseCons(scip, &constraint));

    Edge w_cycle_edge = tp.second;
    SCIP_CALL(SCIPcreateConsBasicLinear(scip, &constraint, "Separation Constraint. W cycle", 1, &w_cycle_edge.var, &coefficient, 0.0, 0.0));
    SCIP_CALL(SCIPaddCons(scip, constraint));
    
    SCIP_CALL(SCIPreleaseCons(scip, &constraint));
  }

  return SCIP_OKAY;
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
    SCIP_VAR **variables = NULL;    

    // Count number of verticies which is at this point equal to 
    // the number of edges
    int number_of_verticies_in_cycle = sb_count(current_cycle);

    // Add doubled edges to the current cycle
    for (int same_edges_index = 0; same_edges_index < sb_count(same_edges); same_edges_index++) {
      Tuple<Edge> edges = same_edges[same_edges_index];
      for (int current_edge_index = 0; current_edge_index < sb_count(current_cycle); current_edge_index++) {
        Edge current_cycle_edge = current_cycle[current_edge_index];
        if (graph_type == TypeOfGraph::DIRECTED) {
          // Branch for directed graphs
          if (are_edges_the_same_directed_graphs(edges.first, current_cycle_edge)) {
            if (current_cycle_edge.var == edges.first.var) {
              sb_push(current_cycle, edges.second);
            } else {
              sb_push(current_cycle, edges.first);
            }
            break;  
          }
        } else {
          // Branch for undirected graphs
          if (are_edges_the_same_undirected_graphs(edges.first, current_cycle_edge)) {
            if (current_cycle_edge.var == edges.first.var) {
              sb_push(current_cycle, edges.second);
            } else {
              sb_push(current_cycle, edges.first);
            }
            break;
          }
        }
      }
    }

    // Add variables to the constraints
    for (int i = 0; i < sb_count(current_cycle); i++) {
      sb_push(variables, current_cycle[i].var);
    }

    // Set values for the variables
    double* values = (double*)allocate_and_zero(sb_count(variables) * sizeof(double));
    for (int i = 0; i < sb_count(variables); i++) {
      values[i] = 1.0;
    }

    // Constraint with less
    SCIP_CALL(SCIPcreateConsBasicLinear(scip, &constraint, "Constraint Less", 
              sb_count(variables), variables, values, -SCIPinfinity(scip), number_of_verticies_in_cycle - 1));
    SCIP_CALL(SCIPaddCons(scip, constraint));

    SCIP_CALL(SCIPreleaseCons(scip, &constraint));

    // Constraint with greater
    SCIP_CALL(SCIPcreateConsBasicLinear(scip, &constraint, "Constraint Greater", 
              sb_count(variables), variables, values, sb_count(current_cycle) - number_of_verticies_in_cycle + 1, SCIPinfinity(scip)));
    SCIP_CALL(SCIPaddCons(scip, constraint));

    SCIP_CALL(SCIPreleaseCons(scip, &constraint));
    
    // Current cycle is Hamilton cycle
    if (number_of_verticies_in_cycle == sb_count(graph)) {
      sb_free(current_cycle);
      sb_free(variables);
      deallocate_and_null((void**)&values);
      return true;
    }

    sb_free(current_cycle);
    sb_free(variables);
    deallocate_and_null((void**)&values);
  }

  return false;
}

/// Analyzes Z-graph and W-graph 
/// @return If Z-graph and W-graph are new Hamilton Decomposition return true otherwise false
bool cycles_are_new_decomposition(SCIP *scip, Edge *z_graph, Edge *w_graph, Tuple<Edge> *same_edges, TypeOfGraph graph_type) {
  bool z_graph_is_complete_cycle = find_cycles_and_add_to_constraints(scip, z_graph, same_edges, graph_type);
  bool w_graph_is_complete_cycle = find_cycles_and_add_to_constraints(scip, w_graph, same_edges, graph_type);
  
  if (z_graph_is_complete_cycle && w_graph_is_complete_cycle) {
    set_text_green();
    printf("Found new Hamilton decomposition\n");
    reset_text_color();
    // print_edges(z_graph);
    // print_edges(w_graph);
    return true;
  } else {
    return false;
  }
}