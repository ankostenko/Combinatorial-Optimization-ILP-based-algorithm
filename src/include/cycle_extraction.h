#pragma once

#include "graph.h"

struct BundleOfGraphs {
  Graph first;
  Graph second;
};

BundleOfGraphs read_cycles_from_stream(FILE* stream, size_t total_size_of_stream);
BundleOfGraphs read_cycles_from_stream(FILE* stream, size_t total_size_of_stream);