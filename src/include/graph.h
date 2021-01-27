#pragma once

#include "stretchy_buffer.h"
#include "scip/scip.h"

struct Vertex {
  int number;
};

struct Edge {
  Vertex start;
  Vertex end;

  SCIP_VAR *var;
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