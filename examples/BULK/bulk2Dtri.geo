// bulk semiconductor in 2D
//
// all measures in um

// in-plane characteristic lengths
la = 0.05;

// device dimension
r = 1;
x = 0.5;
cnt = 0.25;
y = 0.5*r;

// z > 0 part
p1 = newp; Point(p1) = {-x, -y, 0, la};
p2 = newp; Point(p2) = {x, -y, 0, la};
p3 = newp; Point(p3) = {x, y, 0, la};
p4 = newp; Point(p4) = {cnt, y, 0, la / 10};
p5 = newp; Point(p5) = {-cnt, y, 0, la / 10};
p6 = newp; Point(p6) = {-x, y, 0, la};
l1 = newl; Line(l1) = {p1, p2};
l2 = newl; Line(l2) = {p2, p3};
l3 = newl; Line(l3) = {p3, p4};
l4 = newl; Line(l4) = {p4, p5};
l5 = newl; Line(l5) = {p5, p6};
l6 = newl; Line(l6) = {p6, p1};
Line Loop(7) = {6,1,2,3,4,5};
Plane Surface(1) = {7};


Physical Surface(1) = {1};
Physical Line(1) = {4};
Physical Line(2) = {1};
