L = 1;
d = 0.01;
h = 0.1;

Point(1) = {0, 0, 0, d};
Point(2) = {L, 0, 0, d};
Point(3) = {-0.1, 0, 0, d};
Point(4) = {L+0.1, 0, 0, d};

Line(1) = {1, 2};
Line(2) = {3, 1};
Line(3) = {2, 4};

Extrude {0, h, 0} { Line{1,2,3}; Layers {5}; Recombine; }
Extrude {0, 4*h, 0} { Line{8,4,12}; Layers {20}; Recombine; }

Physical Surface("bulk") = {7, 23};
Physical Surface("ndop") = {11, 19};
Physical Surface("pdop") = {15, 27};
Physical Line("anode") = {14};
Physical Line("cathode") = {9,17};
