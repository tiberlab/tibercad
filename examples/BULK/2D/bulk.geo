// bulk semiconductor in 2D
//
// all measures in um

// in-plane characteristic lengths
la = 0.2;
// number of vertical layers
N = 4;

// device dimension
r = 1.3;
x = 0.5;
y = 0.5*r;

// z > 0 part
p1 = newp; Point(p1) = {-x, -y, 0, la};
p2 = newp; Point(p2) = {x, -y, 0, la};
l1 = newl; Line(l1) = {p1, p2};

Extrude Line {1, {0.0, 2*y, 0.0}}
	{Recombine; Layers{N, 1, 1};};

Physical Surface(1) = {1};
Physical Line(1) = {2};
Physical Line(2) = {1};
