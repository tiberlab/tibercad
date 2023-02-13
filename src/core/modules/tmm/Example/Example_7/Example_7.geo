L_AIR = 100;
L_Ag = 20;
L_GaN = 500;
L_Sapphire = 500;

Point(1) = {0,0,0,1};
Point(2) = {L_AIR,0,0,1};
Point(3) = {L_Ag+L_AIR,0,0,1};
Point(4) = {L_GaN/2+L_Ag+L_AIR,0,0,1};
Point(5) = {L_GaN+L_Ag+L_AIR,0,0,1};
Point(6) = {L_Sapphire+L_GaN+L_Ag+L_AIR,0,0,1};
Point(7) = {L_AIR+L_Sapphire+L_GaN+L_Ag+L_AIR,0,0,1};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};
Line(6) = {6,7};

Physical Line("AIR") = {1};
Physical Line("Ag") = {2};
Physical Line("GaN") = {3,4};
Physical Line("Sapphire") = {5,6};

Physical Point("dipole") = {5};









