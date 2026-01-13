L_AIR = 20;
L_MgF2 = 90;
L_Glass = 1e6;
dummy_point = L_Glass/2;
L_ITO = 150;
L_PTAA = 30;
L_PSK = 700;
L_C60 = 30;
L_Ag = 100;



Point(1) = {0,0,0,1};
Point(2) = {L_AIR,0,0,1};
Point(3) = {L_MgF2 + L_AIR,0,0,1};
Point(4) = {dummy_point + L_MgF2 + L_AIR,0,0,1e4};
Point(5) = {L_Glass + L_MgF2 + L_AIR,0,0,1};
Point(6) = {L_ITO + L_Glass + L_MgF2 + L_AIR,0,0,1};
Point(7) = {L_PTAA + L_ITO + L_Glass + L_MgF2 + L_AIR,0,0,1};
Point(8) = {L_PSK + L_PTAA + L_ITO + L_Glass + L_MgF2 + L_AIR,0,0,1};
Point(9) = {L_C60 + L_PSK + L_PTAA + L_ITO + L_Glass + L_MgF2 + L_AIR,0,0,1};
Point(10) = {L_Ag + L_C60 + L_PSK + L_PTAA + L_ITO + L_Glass + L_MgF2 + L_AIR,0,0,1};
Point(11) = {L_AIR + L_Ag + L_C60 + L_PSK + L_PTAA + L_ITO + L_Glass + L_MgF2 + L_AIR,0,0,1};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};
Line(6) = {6,7};
Line(7) = {7,8};
Line(8) = {8,9};
Line(9) = {9,10};
Line(10) = {10,11};

Physical Line("AIR") = {1,10};
Physical Line("MgF2") = {2};
Physical Line("Glass") = {3,4};
Physical Line("ITO") = {5};
Physical Line("PTAA") = {6};
Physical Line("PSK") = {7};
Physical Line("C60") = {8};
Physical Line("Ag") = {9};


Physical Point("incident_wave") = {1};
Physical Point("cathode") = {6};
Physical Point("anode") = {9};



