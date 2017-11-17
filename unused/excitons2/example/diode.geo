Point(1) = {-25,0,0,0.5};
Point(2) = {0,0,0,0.002};
Point(3) = {25,0,0,0.5};
Line(1) = {1,2};
Line(2) = {2,3};


Physical Line("n_side") = {2};
Physical Line("p_side") = {1};
Physical Point("anode") = {1};
Physical Point("cathode") = {3};
