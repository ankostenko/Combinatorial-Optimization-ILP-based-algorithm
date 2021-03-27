# Combinatorial Optimization ILP based algorithm

This algorithm is an implementation of a hamiltonian decomposition algorithm based on an integer linear programming presented in this paper: 

https://arxiv.org/abs/2102.12242

# Problem formulation

INSTANCE. Let *x* and *y* be two Hamiltonian cycles.

QUESTION. Does the multigraph *x* &cup; *y* contain a pair of edge-disjoint Hamiltonian cycles *z* and *w* different from *x* and *y*?

# Results

The algorithm shows comparable results on undirected graphs with known heuristics and significantly better results on directed graphs.
