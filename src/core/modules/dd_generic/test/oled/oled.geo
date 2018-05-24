L1 = 100;
L2 = 100;
d = 1.0;

Point(1) = {0       , 0, 0, d};
Point(2) = {L1      , 0, 0, d};
Point(3) = {L1+L2   , 0, 0, d};
Point(4) = {L1+L2+L1, 0, 0, d};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};

Physical Line("htl") = {1};
Physical Line("eml") = {2};
Physical Line("etl") = {3};
Physical Point("anode") = {1};
Physical Point("cathode") = {4};

