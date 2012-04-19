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

Physical Line (1) = {3};
Physical Line (2) = {1};
Physical Line (3) = {5};
Physical Line (4) = {2};
Physical Line (5) = {4};

Physical Point (1) = {1};
Physical Point (2) = {6};
Physical Point ("int_source") = {2};
Physical Point ("int_drain") = {5};
