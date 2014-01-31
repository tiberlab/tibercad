// Gmsh project created on Thu Feb 21 2012
lc1 = 0.5;
l = 50;
cs = 10.0;
cd = 19.75;
br = 5.0; 


Point(1) = {0, 0, 0, lc1};
Point(2) = {cs, 0, 0, lc1};
Point(3) = {cs + br, 0, 0, lc1};
Point(4) = {cs + br + cd, 0, 0, lc1};
Point(7) = {-l, 0, 0, lc1};
Point(8) = {cs + br + cd + l, 0, 0, lc1};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(6) = {7,1};
Line(7) = {4, 8}; 


Physical Line("source_region") = {1};
Physical Line("barrier1") = {2};
Physical Line("drain_region") = {3};
Physical Line("dd_source_region") = {6};
Physical Line("dd_drain_region") = {7};

Physical Point("source") = {7};
Physical Point("c_source") = {1};
Physical Point("c_drain") = {4};
Physical Point("drain") = {8};

