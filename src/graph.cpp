#include "graph.h"

#include "stretchy_buffer.h"
#include "config_param_parsing.h"

bool are_edges_the_same_undirected_graphs(Edge edge1, Edge edge2) {
  return ((edge1.start.number == edge2.start.number) && (edge1.end.number == edge2.end.number)) || 
         ((edge1.start.number == edge2.end.number) && (edge1.end.number == edge2.start.number));
}

/// Compare two edges for equality
/// @param edge1 first edge
/// @param edge2 second edge
/// @return true if edges are equal false otherwise
bool are_edges_the_same_directed_graphs(Edge edge1, Edge edge2) {
  return (edge1.start.number == edge2.start.number) && (edge1.end.number == edge2.end.number);
}

/// Check if given edges are conjunct for directed graphs
/// @param main_edge Edge that is fixed
/// @param compare_to Edge that is given for comparison
/// @return True if edges are conjunct otherwise return false
bool edges_conjunct_directed_graphs(Edge main_edge, Edge compare_to) {
  return (main_edge.start.number == compare_to.end.number) || 
         (main_edge.end.number == compare_to.start.number);
}

/// Check if given edges are conjunct for undirected graphs
/// @param main_edge Edge that is fixed
/// @param compare_to Edge that is given for comparison
/// @return True if edges are conjunct otherwise return false
bool edges_conjunct_undirected_graphs(Edge main_edge, Edge compare_to) {
  return (main_edge.start.number == compare_to.start.number) ||
         (main_edge.start.number == compare_to.end.number)   || 
         (main_edge.end.number == compare_to.start.number)   ||
         (main_edge.end.number == compare_to.end.number);
}

/// Find cycle
/// @param graph Graph to take cycles from
/// @param graph_type Type of graph directed or undirected
/// @return Newly allocated cycle
Edge* find_cycle(Edge *graph, TypeOfGraph graph_type) {
  Edge *cycle = NULL;
  Edge *current_edge = NULL;

  // Find starting edge  
  for (int i = 0; i < sb_count(graph);i++) {
    current_edge = &graph[i];
    if (current_edge->visited == false) { 
      current_edge->visited = true;
      sb_push(cycle, *current_edge);
      break;
    }
  }

  for (int i = 0; i < sb_count(graph); i++) {
    Edge *next_edge = &graph[i];

    if (next_edge->visited == true) { continue; }
    
    if (graph_type == TypeOfGraph::UNDIRECTED) {
      if (edges_conjunct_undirected_graphs(*current_edge, *next_edge)) {
        // Edges are conjuncted so we add them to the cycle
        sb_push(cycle, *next_edge);
        current_edge = next_edge;
        current_edge->visited = true;
        i = 0;
      }
    } else {
      if (edges_conjunct_directed_graphs(*current_edge, *next_edge)) {
        // Edges are conjuncted so we add them to the cycle
        sb_push(cycle, *next_edge);
        current_edge = next_edge;
        current_edge->visited = true;
        i = 0;
      }
    }
  }

  return cycle;
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

/// Find same edges
/// @param first First cycle
/// @param second Second cycle
/// @param graph_type Type of graph
/// @return Array of same edges for both cycles
Tuple<Edge> *find_same_edges(Graph first, Graph second, TypeOfGraph graph_type) {
  Tuple<Edge> *same_edges = NULL;
  
    for (size_t i = 0; i < first.number_of_edges(); i++) {
      Edge first_cycle_edge = first.edges[i];
      for (size_t j = 0; j < second.number_of_edges(); j++) {
        Edge second_cycle_edge = second.edges[j];

        if (graph_type == TypeOfGraph::DIRECTED) {
          // Check edges for equality
          if (are_edges_the_same_directed_graphs(first_cycle_edge, second_cycle_edge)) {
            sb_push(same_edges, Tuple<Edge>({first_cycle_edge, second_cycle_edge}));
            break;
          }
        } else {
          if (are_edges_the_same_undirected_graphs(first_cycle_edge, second_cycle_edge)) {
            sb_push(same_edges, Tuple<Edge>({first_cycle_edge, second_cycle_edge}));
            break;
          }
        }
      }
    }

  return same_edges;
}