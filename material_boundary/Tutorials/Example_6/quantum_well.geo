Point(1) = {-25, 0, 0, 1};
Point(2) = {-5, 0, 0, 0.1};
Point(3) = {5, 0, 0, 0.1};
Point(4) = {25, 0, 0, 1};

Line (1) = {1, 2};
Line (2) = {2, 3};
Line (3) = {3, 4};


Physical Line  ("well") = {2};
Physical Line  ("buffer1") = {1};
Physical Line  ("buffer2") = {3};


Physical Point ("substrate") = {1};
Physical Point ("top") = {4};
