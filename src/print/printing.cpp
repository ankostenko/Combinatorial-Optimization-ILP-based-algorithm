#include "stdio.h"

#include "printing.h"
#include "graph.h"

void print_cycle_edges(Graph graph) {
  for (size_t i = 0; i < graph.number_of_edges(); i++) {
    printf("%d-%d ", graph.edges[i].start.number, graph.edges[i].end.number);
  }
  printf("\n");
}

void print_edges(Edge *edges) {
  for (int i = 0; i < sb_count(edges); i++) {
    printf("%d %d\n", edges[i].start.number, edges[i].end.number);
  }
  printf("\n");
}

void print_list_of_verticies(std::vector<Vertex> graph) {
  for (Vertex vert : graph) {
    printf("%d ", vert.number);
  }
  printf("\n");
}
