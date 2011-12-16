L = 10;
b=5;
d = 0.1;
d1 = 0.1;

Point(1) = {0, 0, 0, d/10};
Point(2) = {L, 0, 0, d};
Point(3) = {L+b, 0, 0, d1/10};

Line(1) = {1, 2};
Line(2) = {2, 3};

Physical Line("TiO2") = {1};
Physical Line("electrolyte") = {2};
Physical Point("anode") = {1};
Physical Point("cathode") = {3};
