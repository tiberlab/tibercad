//Units: A

// lattice constant
a = 0.3189;
c = 0.5185;

// # repetitions
N = 2.0;

SetFactory("OpenCASCADE");
Box(1) = {0, 0, 0, N*c, N*Sqrt(3)*a, 2*N*a};

Physical Volume("dot") = {1};
