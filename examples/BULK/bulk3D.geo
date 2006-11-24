// bulk semiconductor in 3D
//
// all measures in um

// in-plane characteristic lengths
la = 0.25;
// number of vertical layers
N = 50;

// device dimension
x = 0.5;
y = 0.5;
z = 0.5;

// z > 0 part
p1 = newp; Point(p1) = {-x, -y, -z, la};
p2 = newp; Point(p2) = {-x, y, -z, la};
p3 = newp; Point(p3) = {x, y, -z, la};
p4 = newp; Point(p4) = {x, -y, -z, la};
l1 = newl; Line(l1) = {p1, p4};
l2 = newl; Line(l2) = {p4, p3};
l3 = newl; Line(l3) = {p3, p2};
l4 = newl; Line(l4) = {p2, p1};
Line Loop(1) = {l1, l2, l3, l4};
Plane Surface(1) = {1};
//Recombine Surface {1};

Extrude Surface {1, {0.0, 0.0, 2*z}}
	{Recombine; Layers{N, 1, 1};};

Physical Volume(1) = {1};
Physical Surface(1) = {26};
Physical Surface(2) = {1};
