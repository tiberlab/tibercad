//Units: A

// lattice constant
a = 5.65325;

// # repetitions
N = 2.0;

// should be larger than the piece of material to be constructed
length = N*a;

lc = 0.2;


Point(1) = {-length/2, 0, 0, lc};
Point(2) = { length/2, 0, 0, lc};


Line(1) = {1,2};
Physical Line("dot") = {1};

