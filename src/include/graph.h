#pragma once

#include <vector>

#include "stretchy_buffer.h"
#include "scip/scip.h"

enum TypeOfGraph : int8_t {
  DIRECTED, 
  UNDIRECTED
};

enum GraphName : int8_t {
  Z_GRAPH,
  W_GRAPH,
};

struct Vertex {
  int16_t number;
};

struct Edge {
  SCIP_VAR *var;

  Vertex start;
  Vertex end;
  
  GraphName graph_name;
  bool visited;
  bool visited_for_cycle_search;
  bool fixed;
};

template<typename T> struct Tuple {
  T first;
  T second;
};

struct Graph {
  Vertex *verticies = NULL;
  Edge *edges = NULL;

  size_t number_of_edges() {
    return sb_count(edges);
  }

  void add_edge(Edge edge) {
    sb_push(edges, edge);
  }

  size_t number_of_verticies() {
    return sb_count(verticies);
  }

  void add_vertex(Vertex vert) {
    sb_push(verticies, vert);
  }
};

bool are_edges_the_same_undirected_graphs(Edge edge1, Edge edge2);
bool are_edges_the_same_directed_graphs(Edge edge1, Edge edge2);
bool edges_conjunct_directed_graphs(Edge main_edge, Edge compare_to);
bool edges_conjunct_undirected_graphs(Edge main_edge, Edge compare_to);
Tuple<Edge> *find_same_edges(Graph first, Graph second, TypeOfGraph graph_type);
Edge* find_cycle(Edge *graph, TypeOfGraph graph_type);
bool contains_edge(Edge edge, Edge *edges, TypeOfGraph graph_type);
int degree_of_vertex_in_multigraph(int vertex_index, const std::vector<std::vector<Edge>> &multigraph);
int find_number_of_cycles_in_graph_from_multigraph(std::vector<std::vector<Edge>> &multigraph, GraphName graph_name, TypeOfGraph graph_type);
void visit_edge_in_multigraph(std::vector<std::vector<Edge>> &multigraph, Edge edge, TypeOfGraph graph_type);
int find_number_of_cycles_in_graph(Edge *graph, TypeOfGraph graph_type);
std::tuple<Edge*, Edge*> convert_multigraph_to_two_graph(std::vector<std::vector<Edge>> multigraph, TypeOfGraph graph_type);
void unfix_edges_in_multigraph(std::vector<std::vector<Edge>> &multigraph);
void fix_doubled_edges(std::vector<std::vector<Edge>> &multigraph, Tuple<Edge> *same_edges, TypeOfGraph graph_type);
void fix_edge_in_multigraph(std::vector<std::vector<Edge>> &multigraph, std::vector<int> &broken_verticies, Edge edge, TypeOfGraph graph_type, GraphName graph_name);
void set_brokeness_of_verticies(const std::vector<Edge> &edges, int vertex_number, std::vector<int> &broken_verticies);
bool no_three_incident_fixed_edges(const std::vector<int> &broken_verticies, const std::vector<std::vector<Edge>> &multigraph);