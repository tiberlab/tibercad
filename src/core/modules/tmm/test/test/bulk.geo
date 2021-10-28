Point(1) = {0,0,0,0.2};
Point(2) = {10,0,0,0.2};
Point(3) = {20,0,0,0.2};
Point(4) = {30,0,0,0.2};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};

Physical Line("left") = {1};
Physical Line("middle") = {2};
Physical Line("right") = {3};

Physical Point("point1") = {1};
Physical Point("point2") = {2};
Physical Point("point3") = {3};
Physical Point("point4") = {4};
