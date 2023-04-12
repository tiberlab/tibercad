L_Air = 400;
L_GaN = 1000;


Point(1) = {0,0,0,1};
Point(2) = {L_Air,0,0,1};
Point(3) = {L_Air+L_GaN,0,0,1};

Line(1) = {1,2};
Line(2) = {2,3};

Physical Line("Air") = {1};
Physical Line("GaN") = {2};









