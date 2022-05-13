L_AIR = 20;
L_Glass = 1e6;
dummy_point = L_Glass/2;
L_ITO = 140;
L_PTAA = 20;
L_PSK = 400;
L_SnO2 = 20;
L_Ag = 100;


Point(1) = {0,0,0,1};
Point(2) = {L_AIR,0,0,1};
Point(3) = {dummy_point + L_AIR,0,0,1e4};
Point(4) = {L_Glass + L_AIR,0,0,1};
Point(5) = {L_ITO + L_Glass + L_AIR,0,0,1};
Point(6) = {L_PTAA + L_ITO + L_Glass + L_AIR,0,0,1};
Point(7) = {L_PSK + L_PTAA + L_ITO + L_Glass + L_AIR,0,0,1};
Point(8) = {L_SnO2 + L_PSK + L_PTAA + L_ITO + L_Glass + L_AIR,0,0,1};
Point(9) = {L_Ag + L_SnO2 + L_PSK + L_PTAA + L_ITO + L_Glass + L_AIR,0,0,1};
Point(10) = {L_AIR + L_Ag + L_SnO2 + L_PSK + L_PTAA + L_ITO + L_Glass + L_AIR,0,0,1};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};
Line(6) = {6,7};
Line(7) = {7,8};
Line(8) = {8,9};
Line(9) = {9,10};

Physical Line("AIR") = {1,9};
Physical Line("Glass") = {2,3};
Physical Line("ITO") = {4};
Physical Line("PTAA") = {5};
Physical Line("PSK") = {6};
Physical Line("SnO2") = {7};
Physical Line("Ag") = {8};


Physical Point("point1") = {1};





