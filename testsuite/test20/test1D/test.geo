// Gmsh project created on Thu Feb 21 2012
lc1 = 0.1;
l = 50;
c = 10;
br = 4;
we = 7;


Point(1) = {0, 0, 0, lc1};
Point(2) = {c, 0, 0, lc1};
Point(3) = {c + br, 0, 0, lc1};
Point(4) = {c + br + we, 0, 0, lc1};
Point(5) = {c + br + we + br, 0, 0, lc1};
Point(6) = {c + br + we + br + c, 0, 0, lc1};
Point(7) = {-l, 0, 0, lc1};
Point(8) = {c + br + we + br + c + l, 0, 0, lc1};


Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};
Line(6) = {7,1};
Line(7) = {6,8};


Physical Line("source_region") = {1};
Physical Line("barrier1") = {2};
Physical Line("well") = {3};
Physical Line("barrier2") = {4};
Physical Line("drain_region") = {5};
Physical Line("dd_source_region") = {6};
Physical Line("dd_drain_region") = {7};

Physical Point("source") = {7};
Physical Point("int_source") = {1};
Physical Point("int_drain") = {6};
Physical Point("drain") = {8};

