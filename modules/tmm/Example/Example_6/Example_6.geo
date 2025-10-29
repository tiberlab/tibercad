L_AIR1 = 100;
L_MgF2 = 90;
L_Glass = 1e6;
dummy_point = L_Glass/2;
L_ITO = 130;
L_PTAA = 10;
L_PSK = 700;
L_C60 = 30;
L_Ag = 100;
L_AIR2 = 100;

x0 = -(L_Glass + L_MgF2 + L_AIR1);

Point(1) = {x0,0,0,1};
Point(2) = {x0+L_AIR1,0,0,1};
Point(3) = {x0+L_MgF2 + L_AIR1,0,0,1};
Point(4) = {x0+dummy_point + L_MgF2 + L_AIR1,0,0,L_Glass/100};
Point(5) = {x0+L_Glass + L_MgF2 + L_AIR1,0,0,1};
Point(6) = {x0+L_ITO + L_Glass + L_MgF2 + L_AIR1,0,0,1};
Point(7) = {x0+L_PTAA + L_ITO + L_Glass + L_MgF2 + L_AIR1,0,0,1};
Point(8) = {x0+L_PSK + L_PTAA + L_ITO + L_Glass + L_MgF2 + L_AIR1,0,0,1};
Point(9) = {x0+L_C60 + L_PSK + L_PTAA + L_ITO + L_Glass + L_MgF2 + L_AIR1,0,0,1};
Point(10) = {x0+L_Ag + L_C60 + L_PSK + L_PTAA + L_ITO + L_Glass + L_MgF2 + L_AIR1,0,0,1};
Point(11) = {x0+L_AIR2 + L_Ag + L_C60 + L_PSK + L_PTAA + L_ITO + L_Glass + L_MgF2 + L_AIR1,0,0,1};

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

Physical Line("AIR1") = {1};
Physical Line("MgF2") = {2};
Physical Line("Glass") = {3,4};
Physical Line("ITO") = {5};
Physical Line("PTAA") = {6};
Physical Line("PSK") = {7};
Physical Line("C60") = {8};
Physical Line("Ag") = {9};
Physical Line("AIR2") = {10};

Physical Point("incident_wave") = {1};
Physical Point("cathode") = {6};
Physical Point("anode") = {9};



