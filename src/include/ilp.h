#pragma once

#include "scip/scip.h"
#include "graph.h"

SCIP_RETCODE add_first_constraint(SCIP *scip, Graph first, Graph second);
SCIP_RETCODE add_second_constraint(SCIP *scip, Graph first, Graph second, TypeOfGraph graph_type);
SCIP_RETCODE add_third_or_fourth_constraint(SCIP *scip, Graph graph, Tuple<Edge> *same_edges, TypeOfGraph graph_type, const char *name_of_constraint);
SCIP_RETCODE add_seventh_constraint(SCIP *scip, Tuple<Edge> *same_edges);
bool cycles_are_new_decomposition(SCIP *scip, Edge *z_graph, Edge *w_graph, Tuple<Edge> *same_edges, TypeOfGraph graph_type);