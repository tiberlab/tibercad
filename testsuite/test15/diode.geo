scale = 20;

Point(1) = {-50 * scale,0,0,0.5 * scale};
Point(2) = {0,0,0,0.002 * scale};
Point(3) = {50 *scale ,0,0,0.5 * scale};
Line(1) = {1,2};
Line(2) = {2,3};



Physical Line("p_side") = {1};
Physical Line("n_side") = {2};
Physical Point("anode") = {1};
Physical Point("cathode") = {3};

