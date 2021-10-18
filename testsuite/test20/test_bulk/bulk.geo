L = 10;
d = 0.1;

Point(1) = {0, 0, 0, d};
Point(2) = {L, 0, 0, d};

Line(1) = {1, 2};

Physical Line("bulk") = {1};
Physical Point("anode") = {1};
Physical Point("cathode") = {2};
