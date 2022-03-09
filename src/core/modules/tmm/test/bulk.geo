L_AIR1 = 20;
L_Glass = 100;
L_ITO = 142;
L_PTAA = 20;
L_PSK = 400;
L_SnO2 = 20;
L_Ag = 100;
L_AIR2 = 10;

Point(1) = {0,0,0,1};
Point(2) = {L_AIR1,0,0,1};
Point(3) = {L_Glass + L_AIR1,0,0,1};
Point(4) = {L_ITO + L_Glass + L_AIR1,0,0,1};
Point(5) = {L_PTAA + L_ITO + L_Glass + L_AIR1,0,0,1};
Point(6) = {L_PSK + L_PTAA + L_ITO + L_Glass + L_AIR1,0,0,1};
Point(7) = {L_SnO2 + L_PSK + L_PTAA + L_ITO + L_Glass + L_AIR1,0,0,1};
Point(8) = {L_Ag + L_SnO2 + L_PSK + L_PTAA + L_ITO + L_Glass + L_AIR1,0,0,1};
Point(9) = {L_AIR2 + L_Ag + L_SnO2 + L_PSK + L_PTAA + L_ITO + L_Glass + L_AIR1,0,0,1};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};
Line(6) = {6,7};
Line(7) = {7,8};
Line(8) = {8,9};

Physical Line("AIR1") = {1};
Physical Line("Glass") = {2};
Physical Line("ITO") = {3};
Physical Line("PTAA") = {4};
Physical Line("PSK") = {5};
Physical Line("SnO2") = {6};
Physical Line("Ag") = {7};
Physical Line("AIR2") = {8};

Physical Point("point1") = {1};
Physical Point("point2") = {2};
Physical Point("point3") = {3};
Physical Point("point4") = {4};
Physical Point("point5") = {5};
Physical Point("point6") = {6};
Physical Point("point7") = {7};
Physical Point("point8") = {8};
Physical Point("point9") = {9};




