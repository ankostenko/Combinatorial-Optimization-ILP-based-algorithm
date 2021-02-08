#include <string>
#include <sstream>

#include "stretchy_buffer.h"
#include "cycle_extraction.h"
#include "utils.h"

/// Extracts verticies from a string
/// @param str string to extract verticies from
/// @return Graph cycle of parsed verticies
Graph extract_verticies_from_string(const char *str) {
  Graph cycle;

  std::stringstream ss(str);
	std::string token;
  while (std::getline(ss, token, ' ')) {
    if (token == "\n") { continue; }
    cycle.add_vertex(Vertex{ atoi(token.c_str()) });
  }

  return cycle;
}

/// Create edges from extracted verticies
/// @param graph To take verticies from
/// @return Return same graph with newly constructed edges
Graph create_edges_from_verticies(Graph *graph) {
  for (size_t i = 0; i < graph->number_of_verticies(); i++) {
    // If we iterated to the last vertex
    if (i + 1 == graph->number_of_verticies()) {
      graph->add_edge(Edge{ .start = graph->verticies[i], .end = graph->verticies[0] });
    } else {
      graph->add_edge(Edge{ .start = graph->verticies[i], .end = graph->verticies[i + 1] });
    }
  }
  return *graph;
}

/// Read two cycles from the stream and parse it to Graph structures
/// @param bundle bundle to fill in
/// @param stream stream to read from
/// @param total_size_of_stream maximum size of a file or any other stream
/// @return True if read is successful otherwise return false
bool read_cycles_from_stream(BundleOfGraphs *bundle, FILE* stream, size_t total_size_of_stream) {
  Graph first = { };
  Graph second = { };

  // Read first cycle
  char *first_cycle_line = (char*)allocate_and_zero(total_size_of_stream);
  if (fgets(first_cycle_line, total_size_of_stream, stream) == NULL) {
    return false;
  }
  first = extract_verticies_from_string(first_cycle_line);
  first = create_edges_from_verticies(&first);

  // Read second cycle
  char *second_cycle_line = (char*)allocate_and_zero(total_size_of_stream);
  if (fgets(second_cycle_line, total_size_of_stream, stream) == NULL) {
    report_problem_to_user(true, "Unexpected end of file there's no second cycle");
    return false;
  }
  second = extract_verticies_from_string(second_cycle_line);
  second = create_edges_from_verticies(&second);

  deallocate_and_null((void**)&first_cycle_line);
  deallocate_and_null((void**)&second_cycle_line);

  *bundle = { first, second };

  // Skip new line
  if (fgetc(stream) == EOF) {
    return false;
  }

  return true;
}