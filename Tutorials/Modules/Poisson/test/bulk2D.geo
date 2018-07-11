L = 2;
H = 1;
d = 0.03;
dc = 0.01;

Point(1) = {0, 0, 0, d};
Point(2) = {L, 0, 0, d};
Point(3) = {L, H/2, 0, d};
Point(4) = {3*L/4, H/2, 0, dc}; 
Point(5) = {L/2, H/2, 0, dc}; 
Point(6) = {L/2, H, 0, d};
Point(7) = {L/4, H, 0, dc};
Point(8) = {0, H, 0, d}; 

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 1};

Line Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8};
Plane Surface(1) = {1};

Physical Surface("bulk") = {1};
Physical Line("anode") = {6, 7};
Physical Line("cathode") = {2};
