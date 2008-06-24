// bulk semiconductor in 3D
//
// all measures in um

// in-plane characteristic lengths
la = 0.25;
lct = 0.025;
// number of vertical layers
N = 20;

// device dimension (half of side lengths)
x = 0.5;
y = 0.5;
// height
z = 0.5;

// the same for the contact
x_ct = 0.1;
y_ct = 0.1;

// z > 0 part
p1 = newp; Point(p1) = {-x, -y, 0, la};
p2 = newp; Point(p2) = {-x, y, 0, la};
p3 = newp; Point(p3) = {x, y, 0, la};
p4 = newp; Point(p4) = {x, -y, 0, la};

p5 = newp; Point(p5) = {-x_ct, -y_ct, 0, lct};
p6 = newp; Point(p6) = {-x_ct, y_ct, 0, lct};
p7 = newp; Point(p7) = {x_ct, y_ct, 0, lct};
p8 = newp; Point(p8) = {x_ct, -y_ct, 0, lct};

l1 = newl; Line(l1) = {p1, p4};
l2 = newl; Line(l2) = {p4, p3};
l3 = newl; Line(l3) = {p3, p2};
l4 = newl; Line(l4) = {p2, p1};

l5 = newl; Line(l5) = {p5, p8};
l6 = newl; Line(l6) = {p8, p7};
l7 = newl; Line(l7) = {p7, p6};
l8 = newl; Line(l8) = {p6, p5};

Line Loop(1) = {l1, l2, l3, l4};
Line Loop(2) = {l5, l6, l7, l8};
Plane Surface(1) = {1,2};
Plane Surface(2) = {2};

//Recombine Surface {1};
//Recombine Surface {2};

Extrude Surface {1, {0.0, 0.0, z}}
	{Recombine; Layers{N, 1, 1};};
Extrude Surface {2, {0.0, 0.0, z}}
	{Recombine; Layers{N, 1, 1};};

anode = news - 1;
Physical Volume(1) = {1,2};
Physical Surface(1) = {anode};
Physical Surface(2) = {1, 2};
