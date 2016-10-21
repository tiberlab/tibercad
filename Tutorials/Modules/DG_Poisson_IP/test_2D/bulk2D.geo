L = 2;
H = 1;
d = 0.01;

Point(1) = {0, 0, 0, d};
Point(2) = {L, 0, 0, d};
Point(3) = {L, H, 0, d};
Point(4) = {0, H, 0, d}; 

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

Physical Surface("bulk") = {1};
Physical Line("anode") = {3};
Physical Line("cathode") = {1};

