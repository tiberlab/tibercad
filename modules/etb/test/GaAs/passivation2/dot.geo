//Units: A

// lattice constant
a = 0.565325;

// # repetitions
N = 2.0;

SetFactory("OpenCASCADE");
Box(1) = {0, 0, 0, N*a, N*a, N*a};

Physical Volume("dot") = {1};
