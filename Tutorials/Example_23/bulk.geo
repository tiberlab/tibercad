L1 = 2;
L2 = 196;
L3 = 2;
d1 = 0.5;
d2 = 0.1;

Point(1) = {0       , 0, 0, d2};
Point(2) = {L1      , 0, 0, d1};
Point(3) = {L1+L2   , 0, 0, d1};
Point(4) = {L1+L2+L3, 0, 0, d2};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};

Physical Line ("bulk"       ) = {2};
Physical Line ("int_anode"  ) = {1};
Physical Line ("int_cathode") = {3};
Physical Point("anode"      ) = {1};
Physical Point("cathode"    ) = {4};
