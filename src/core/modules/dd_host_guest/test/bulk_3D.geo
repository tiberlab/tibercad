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

Extrude {0, h, 0} { Line{1,2,3}; Layers {5};  }
Extrude {0, 4*h, 0} { Line{8,4,12}; Layers {20};  }

Extrude {0, 0, h} { Surface{7, 23, 11, 19, 15, 27}; Layers {5};  }
Extrude {0, 0, 4*h} { Surface{71, 49, 159, 115, 93, 137}; Layers{20};  }

Physical Volume("bulk") = {2, 7, 1, 8};
Physical Volume("ndop") = {4, 11, 10, 3};
Physical Volume("pdop") = {5, 6, 9, 12};
Physical Surface("anode") = {282, 128, 150};
Physical Surface("cathode") = {114, 246, 268, 92};
