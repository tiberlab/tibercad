cl = 2;
Point(1) = {0,0,0,cl};
Point(2) = {342,0,0,cl};
Point(3) = {342,156,0,cl};
Point(4) = {0,156,0,cl};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Line Loop(1) = {1,2,3,4};
Plane Surface(1) = {1};

Physical Surface("bulk") = {1};
