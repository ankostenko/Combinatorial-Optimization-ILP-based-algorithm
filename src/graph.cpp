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

/// Check if edge contains in 'edges' array
bool contains_edge(Edge edge, Edge *edges, TypeOfGraph graph_type) {
  for (int i = 0; i < sb_count(edges); i++) {
    if (TypeOfGraph::DIRECTED == graph_type) {
      if (are_edges_the_same_directed_graphs(edges[i], edge)) {
        return true;
      }
    } else {
      if (are_edges_the_same_undirected_graphs(edges[i], edge)) {
        return true;
      }
    }
  }

  return false;
}

/// Calculate degree of vertex of graph Z
/// @param vertex_index Index of vertex
/// @param multigraph List of adjacency of both both W and Z
/// @return Degree of a given vertex
int degree_of_vertex_in_multigraph(int vertex_index, std::vector<std::vector<Edge>> &multigraph) {
  std::vector<Edge> &vec = multigraph[vertex_index - 1];
  int degree_of_vertex = 0;
  for (Edge e: vec) {
    if (e.graph_name == GraphName::Z_GRAPH) {
      degree_of_vertex++;
    }
  }

  return degree_of_vertex;
}

/// Fix edge in multigraph and set in which graph it goes
void visit_edge_in_multigraph(std::vector<std::vector<Edge>> &multigraph, Edge edge, TypeOfGraph graph_type) {
  std::vector<Edge> &start_vector = multigraph[edge.start.number - 1];
  std::vector<Edge> &end_vector = multigraph[edge.end.number - 1];

  for (Edge &e : start_vector) {
    if (TypeOfGraph::DIRECTED == graph_type) {
      if (are_edges_the_same_directed_graphs(e, edge)) {
        e.visited = true;
      }
    } else {
      if (are_edges_the_same_undirected_graphs(e, edge)) {
        e.visited = true;
      }
    }
  }

  for (Edge &e: end_vector) {
    if (TypeOfGraph::DIRECTED == graph_type) {
      if (are_edges_the_same_directed_graphs(e, edge)) {
        e.visited = true;
      }
    } else {
      if (are_edges_the_same_undirected_graphs(e, edge)) {
        e.visited = true;
      }
    }
  }
}


/// Find number of cycles in given graph
int find_number_of_cycles_in_graph_from_multigraph(std::vector<std::vector<Edge>> multigraph, GraphName graph_name, TypeOfGraph graph_type) {
  int number_of_cycles = 0;
  Edge current_vertex;
  size_t i = 0;

  while (true) {
    // Find starting vertex
    for (i = 0; i < multigraph.size(); i++){
      std::vector<Edge> &start_vertex_array = multigraph[0];
      Edge start_vertex = { .start = { -1 } };
      for (Edge &e: start_vertex_array) {
        if (e.graph_name == graph_name && e.visited == false) {
          current_vertex = e;
          start_vertex = e;
          visit_edge_in_multigraph(multigraph, e, graph_type);
          break;
        }
      }
      if (start_vertex.start.number == -1) { 
        return number_of_cycles; 
      } else {
        number_of_cycles += 1;
        break;
      }
    }
    

    while (true) {
      std::vector<Edge> &vec = multigraph[current_vertex.end.number - 1];

      for (i = 0; i < vec.size(); i++) {
        Edge &e = vec[i];
        if (e.visited == true) { continue; }
        if ((e.graph_name == graph_name) && ((e.start.number == current_vertex.end.number) || (e.end.number == current_vertex.end.number))) {
          visit_edge_in_multigraph(multigraph, e, graph_type);
          // Swap start and end of the edge 
          // It should be done because if the previous edge was in another graph
          // Orientation of edge could be different
          if (e.end.number == current_vertex.end.number) {
            std::swap(e.start.number, e.end.number);
          }
          current_vertex = e;
          break;
        }
      }
      if (i == 4) { break; }
    }
  }
}