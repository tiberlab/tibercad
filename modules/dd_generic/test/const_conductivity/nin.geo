L = 10;
d = 0.002;
d2 = 0.001;

Point(1) = {0, 0, 0, d};
Point(2) = {L, 0, 0, d};
Point(3) = {-0.2, 0, 0, d2};
Point(4) = {L+0.2, 0, 0, d2};

Line(1) = {1, 2};
Line(2) = {3, 1};
Line(3) = {2, 4};

Physical Line("intrinsic") = {1};
Physical Line("ndop") = {2};
Physical Line("pdop") = {3};
Physical Point("anode") = {4};
Physical Point("cathode") = {3};
