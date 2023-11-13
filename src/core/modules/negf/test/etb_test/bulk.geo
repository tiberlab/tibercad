// Gmsh project created on Thu Apr 21 2022
lc1 = 0.565;
l = 5*lc1;
cs = 20*lc1;


Point(1) = {0, 0, 0, lc1};
Point(2) = {cs, 0, 0, lc1};
Point(7) = {-l, 0, 0, lc1};
Point(8) = {cs + l, 0, 0, lc1};

Line(1) = {1,2};
Line(2) = {7,1};
Line(3) = {2, 8}; 


Physical Line("bulk") = {1};
Physical Line("dd_source_region") = {2};
Physical Line("dd_drain_region") = {3};

Physical Point("source") = {7};
Physical Point("c_source") = {1};
Physical Point("c_drain") = {2};
Physical Point("drain") = {8};

