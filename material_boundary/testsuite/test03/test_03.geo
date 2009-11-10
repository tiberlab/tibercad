scale = 20;

Point(1) = {-50 * scale,0,0,0.5 * scale};
Point(2) = {-25 *scale,0,0,0.5 * scale};            
Point(3) = {0,0,0,0.002 * scale};
Point(4) = {25 * scale,0,0,0.05 * scale};
Point(5) = {50 *scale ,0,0,0.5 * scale};
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};



Physical Line(1) = {1};
Physical Line(2) = {2};
Physical Line(3) = {3};
Physical Line(4) = {4};



Physical Point(1) = {1};
Physical Point(2) = {2};
Physical Point(3) = {4};
Physical Point(4) = {5};


