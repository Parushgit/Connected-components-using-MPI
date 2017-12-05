# Connected-components-using-MPI

Your task in this assignment is to implement efficient C++/MPI function connected_components. We make several assumptions:

Undirected graph on which we are operating is too large to be represented in the memory of a single compute node.
We have p = q * q ranks available.
The graph has n nodes, and we have that q divides n.
The graph is represented by the adjacency matrix A, in which 1 indicates edge and 0 means no edge.
Adjacency matrix A is 2D-decomposed using q by q row-wise grid of ranks.
Instructions
Download A1.tar.bz2 that contains the backbone of the project (link).
Edit a1.hpp and implement connected_components. You may include extra files and implement additional functionality as needed, however you are not allowed to change the signature of the connected_components. Moreover, you should not edit a1.cpp as this is example of the engine that will be testing your code. Arguments of the connected_components are as follows:
A adjacency matrix, row-wise block of size n/q by n/q.
n total number of nodes in the graph.
q dimension of the rank grid (p = q * q).
out path to the output file where the assignment of nodes to the connected components should be stored (see below).
comm communicator with p = q * q ranks to work with.
Function must return the total number of connected components found.
Output file should be written as binary. Output should store P: vector of ints, where size of P is n, and P[i] is an identifier of the connected component to which node i belongs.
When implementing connected_components you may assume that all arguments are correct, and there is no need to test that. Finally, the project backbone provides a simple Makefile. You can edit it as needed for your project. However, when invoked without arguments, Makefile must produce executable a1 from a1.cpp.
