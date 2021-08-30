//Units: A

// lattice constant
a = 5.65325;

// # repetitions
N = 2.0;

// should be larger than the piece of material to be constructed
well_length = N*a;

lc = 0.2;


Point(1) = {-well_length/2,0,0,lc};
Point(2) = {well_length/2, 0, 0, lc};

//Physical Point("cathode") = {1};
//Physical Point("anode") = {2};

Line(1) = {1,2};
Physical Line("bulk") = {1};

