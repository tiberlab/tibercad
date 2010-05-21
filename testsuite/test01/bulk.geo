L = 1;
d = 0.01;

Point(1) = {0, 0, 0, d};
Point(2) = {L, 0, 0, d};
Point(3) = {-0.1, 0, 0, d};
Point(4) = {L+0.1, 0, 0, d};

Line(1) = {1, 2};
Line(2) = {3, 1};
Line(3) = {2, 4};

Physical Line("bulk") = {1};
Physical Line("ndop") = {2,3};
Physical Point("anode") = {3};
Physical Point("cathode") = {4};
