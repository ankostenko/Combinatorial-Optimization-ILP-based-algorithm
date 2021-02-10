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
int degree_of_vertex_in_multigraph(int vertex_index, const std::vector<std::vector<Edge>> &multigraph) {
  const std::vector<Edge> &vec = multigraph[vertex_index - 1];
  int degree_of_vertex = 0;
  for (const Edge &e: vec) {
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

/// Find number of components in a given graph
/// @param graph Original graph
/// @param graph_type Type of graph
/// @return Number of components in a given graph
int find_number_of_cycles_in_graph(Edge *graph, TypeOfGraph graph_type) {
  int number_of_components_in_graph = 0;
  while (true) {
    Edge *cycle = find_cycle(graph, graph_type);
    if (cycle == NULL) { break; }
    number_of_components_in_graph += 1;
  }

  return number_of_components_in_graph;
}

/// Extract Z graph and W graph from multigraph
/// @param multigraph Original multigraph
/// @param graph_type Type of graph
/// @return Tuple of newly allocated Z and W graphs
std::tuple<Edge*, Edge*> convert_multigraph_to_two_graph(std::vector<std::vector<Edge>> multigraph, TypeOfGraph graph_type) {
  Edge *z_graph = NULL;
  Edge *w_graph = NULL;

  // Fill Z graph and W graph with better solution
  for (std::vector<Edge> vec: multigraph) {
    for (Edge e: vec) {
      // Filling Z graph
      if (e.graph_name == GraphName::Z_GRAPH) {
        if (contains_edge(e, z_graph, graph_type) == false) {
          sb_push(z_graph, e);                
        }
      }

      // Fillign W graph
      if (e.graph_name == GraphName::W_GRAPH) {
        if (contains_edge(e, w_graph, graph_type) == false) {
          sb_push(w_graph, e);
        }
      }
    }
  }

  return std::make_tuple(z_graph, w_graph);
}

/// Check if there's no three incident fixed edges
/// Three incident fixed edges mean that attempt is no longer valid
/// And we need to start again.
/// @param broken_verticies Array of broken verticies
/// @param multigraph Multigraph of all edges
/// @return True if no incident fixed edges otherwise false
bool no_three_incident_fixed_edges(const std::vector<int> &broken_verticies, const std::vector<std::vector<Edge>> &multigraph) {
  // Check if there's no verticies with 3 incident fixed edges
  for (int vertex_index: broken_verticies) {
    const std::vector<Edge> &vec = multigraph[vertex_index - 1];
    int number_of_fixed_edges = 0;
    // Z graph
    for (const Edge &e: vec) {
      if (e.graph_name == GraphName::Z_GRAPH) {
        number_of_fixed_edges += e.fixed;
      }
    }

    // Number of fixed incident edges is three or more so we should stop attempt
    if (number_of_fixed_edges >= 3) { return false; }

    // W graph
    number_of_fixed_edges = 0;
    for (const Edge &e: vec) {
      if (e.graph_name == GraphName::W_GRAPH) {
        number_of_fixed_edges += e.fixed;
      }
    }

    // Number of fixed incident edges is three or more so we should stop attempt
    if (number_of_fixed_edges >= 3) { return false; }
  }
  
  return true;
}

/// Unfix doubled edges in multigraph
/// @param multigraph Original multigraph
/// @param same_edges Array of doubled edges
/// @param graph_type Type of graph
void unfix_edges_in_multigraph(std::vector<std::vector<Edge>> &multigraph) {
  for (auto &vec: multigraph) {
    for (Edge &e: vec) {
      e.fixed = false;
    }
  }
}

/// Fix doubled edges in multigraph
/// @param multigraph Original multigraph
/// @param same_edges Array of doubled edges
/// @param graph_type Type of graph
void fix_doubled_edges(std::vector<std::vector<Edge>> &multigraph, Tuple<Edge> *same_edges, TypeOfGraph graph_type) {
  for (auto &vec: multigraph) {
    for (Edge &e: vec) {
      if (is_edge_contains_in(same_edges, e, graph_type)) {
        e.fixed = true;
      }
    }
  }
}

/// Add or delete vertex from broken verticies array
/// @param edges Which adjacent edges to look in
/// @param vertex_number Number of a given vertex
/// @param broken_verticies Array of broken verticies
void set_brokeness_of_verticies(const std::vector<Edge> &edges, int vertex_number, std::vector<int> &broken_verticies) {
  // Check if the verticies now broken
  // 1. Find degree of the vertex
  int degree_of_vertex = 0;
  for (const Edge &e: edges) {
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
/// @param multigraph Original multigraph
/// @param broken_verticies Array of broken verticies
/// @param edge Edge to fix in multigraph
/// @param graph_type Type of graph
/// @param graph_name Name of graph to fix edge in
void fix_edge_in_multigraph(std::vector<std::vector<Edge>> &multigraph, std::vector<int> &broken_verticies, Edge edge, 
                                                                            TypeOfGraph graph_type, GraphName graph_name) {
  std::vector<Edge> &start_vector = multigraph[edge.start.number - 1];
  std::vector<Edge> &end_vector = multigraph[edge.end.number - 1];

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