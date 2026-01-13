L_Air = 200;
L_GaN = 500;
L_Ag  = 200;

Point(1) = {0,0,0,1};
Point(2) = {L_Air,0,0,1};
Point(3) = {L_Air+L_GaN,0,0,1};
Point(4) = {L_Air+L_GaN+L_Ag,0,0,1};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};

Physical Line("Air") = {1};
Physical Line("GaN") = {2};
Physical Line("Ag") = {3};









