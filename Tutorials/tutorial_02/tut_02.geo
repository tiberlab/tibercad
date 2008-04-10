scale = 20;

Point(1) = {-50 * scale,0,0,0.5 * scale};
Point(2) = {0,0,0,0.002 * scale};
Point(3) = {50 *scale ,0,0,0.5 * scale};
Line(1) = {1,2};
Line(2) = {2,3};



Physical Line(1) = {1};
Physical Line(2) = {2};
Physical Point(1) = {1};
Physical Point(2) = {3};


