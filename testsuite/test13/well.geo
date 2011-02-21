Point(1) = {-100, 0, 0, 4};
Point(2) = {-20, 0, 0, 1.0};
Point(3) = {-10, 0, 0, 0.1};
Point(4) = {10, 0, 0, 0.1};
Point(5) = {20, 0, 0, 1.0};
Point(6) = {100, 0, 0, 4};



Line (1) = {1, 2};
Line (2) = {2, 3};
Line (3) = {3, 4};
Line (4) = {4, 5};
Line (5) = {5, 6};

Physical Point ("left") = {1};
Physical Line ("well") = {1};
Physical Line ("buffer") = {2,3,4,5};
Physical Point ("right") = {6};
