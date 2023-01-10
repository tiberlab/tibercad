L_AIR = 5e-3;
L_Ag = 20e-3;



Point(1) = {0,0,0,1e-3};
Point(2) = {L_AIR,0,0,1e-3};
Point(3) = {L_Ag + L_AIR,0,0,1e-3};
Point(4) = {L_AIR + L_Ag + L_AIR,0,0,1e-3};


Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};

Physical Line("AIR") = {1,3};
Physical Line("Ag") = {2};

Physical Point("point1") = {1};








