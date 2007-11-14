// bulk semiconductor in 2D
//
// all measures in um

// in-plane characteristic lengths
la = 0.1;
// number of vertical layers
N = 20;

// device dimension
x = 2;
cnt = 0.25;
y = 1;

// z > 0 part
p1 = newp; Point(p1) = {-x, 0, 0, la};
p2 = newp; Point(p2) = {-cnt, 0, 0, la / 5};
p3 = newp; Point(p3) = {0, 0, 0, la};
p4 = newp; Point(p4) = {cnt, 0, 0, la / 5};
p5 = newp; Point(p5) = {x, 0, 0, la};
l1 = newl; Line(l1) = {p1, p2};
l2 = newl; Line(l2) = {p2, p3};
l3 = newl; Line(l3) = {p3, p4};
l4 = newl; Line(l4) = {p4, p5};

Extrude Line {l1, {0.0, y, 0.0}}
//	{Recombine; Layers{N, 1, 1};};
	{Layers{N, 1, 1};};
Extrude Line {l2, {0.0, y, 0.0}}
//	{Recombine; Layers{N, 1, 1};};
	{Layers{N, 1, 1};};
Extrude Line {l3, {0.0, y, 0.0}}
//	{Recombine; Layers{N, 1, 1};};
	{Layers{N, 1, 1};};
Extrude Line {l4, {0.0, y, 0.0}}
//	{Recombine; Layers{N, 1, 1};};
	{Layers{N, 1, 1};};

Physical Surface(1) = {1};
Physical Line(1) = {9,13};
Physical Line(2) = {1,2,3,4};
Physical Line(3) = {5,17};
