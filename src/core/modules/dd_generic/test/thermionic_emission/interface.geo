L = 0.1;
d = 0.001;
d2 = 0.01;

Point(1) = {0, 0, 0, d2};
Point(2) = {L, 0, 0, d};
Point(3) = {2*L, 0, 0, d2};

Line(1) = {1, 2};
Line(2) = {2, 3};

Physical Line("left") = {1};
Physical Line("right") = {2};
Physical Point("anode") = {1};
Physical Point("cathode") = {3};
