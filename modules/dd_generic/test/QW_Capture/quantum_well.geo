buf = 15;
well = 5;

Point(1) = {-200, 0, 0, 5};
Point(2) = {-20, 0, 0, 1.0};
Point(3) = {-buf, 0, 0, 1.0};
Point(4) = {-well, 0, 0, 0.1};
Point(5) = {well, 0, 0, 0.1};
Point(6) = {buf, 0, 0, 1.0};
Point(7) = {20, 0, 0, 1.0};
Point(8) = {200, 0, 0, 5};



Line (1) = {1, 2};
Line (2) = {2, 3};
Line (3) = {3, 4};
Line (4) = {4, 5};
Line (5) = {5, 6};
Line (6) = {6, 7};
Line (7) = {7, 8};

Physical Point ("cathode") = {1};

Physical Line ("nside") = {1};
Physical Line ("buffer_classical") = {2,6};
Physical Line ("buffer_quantum") = {3,5};
Physical Line ("well") = {4};
Physical Line ("pside") = {7};

Physical Point ("anode") = {8};
