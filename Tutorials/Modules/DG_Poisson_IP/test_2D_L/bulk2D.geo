L = 2;
H = 1;
d = 0.05;
dc = 0.001;

Point(1) = {0, 0, 0, d};
Point(2) = {L, 0, 0, d};
Point(3) = {L, H/2, 0, d};
Point(4) = {L/2, H/2, 0, dc}; 
Point(5) = {L/2, H, 0, d}; 
Point(6) = {0, H, 0, d}; 

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 1};

Line Loop(1) = {1, 2, 3, 4, 5, 6};
Plane Surface(1) = {1};

Physical Surface("bulk") = {1};
Physical Line("anode") = {5};
Physical Line("cathode") = {1};
