#pragma once

#include "stretchy_buffer.h"
#include "scip/scip.h"

enum TypeOfGraph {
  DIRECTED, 
  UNDIRECTED
};

struct Vertex {
  int number;
};

struct Edge {
  Vertex start;
  Vertex end;

  SCIP_VAR *var;
  bool visited;
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