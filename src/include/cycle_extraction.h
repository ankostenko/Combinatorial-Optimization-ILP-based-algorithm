#pragma once

#include "graph.h"

struct BundleOfGraphs {
  Graph first;
  Graph second;
};

bool read_cycles_from_stream(BundleOfGraphs *bundle, FILE* stream, size_t total_size_of_stream);