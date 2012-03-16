// Gmsh project created on Thu Feb 21 2012
lc1 =0.25;

c = 10;
br = 2;
we = 5;


Point(1) = {0, 0, 0, lc1};
Point(2) = {c, 0, 0, lc1};
Point(3) = {c + br, 0, 0, lc1};
Point(4) = {c + br + we, 0, 0, lc1};
Point(5) = {c + br + we + br, 0, 0, lc1};
Point(6) = {c + br + we + br + c, 0, 0, lc1};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};

Physical Line("source_region") = {1};
Physical Line("barrier1") = {2};
Physical Line("well") = {3};
Physical Line("barrier2") = {4};
Physical Line("drain_region") = {5};

Physical Point("source") = {1};
Physical Point("int_source") = {2};
Physical Point("int_drain") = {5};
Physical Point("drain") = {6};

