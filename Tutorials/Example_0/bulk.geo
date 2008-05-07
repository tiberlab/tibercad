L = 1;
d = 0.01;

Point(1) = {0, 0, 0, d};
Point(2) = {L, 0, 0, d};

Line(1) = {1, 2};

Physical Line(1) = {1};
Physical Point(1) = {1};
Physical Point(2) = {2};
