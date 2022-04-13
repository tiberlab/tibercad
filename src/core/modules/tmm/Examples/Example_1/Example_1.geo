L_AIR = 5;
L_Ag = 20;



Point(1) = {0,0,0,1};
Point(2) = {L_AIR,0,0,1};
Point(3) = {L_Ag + L_AIR,0,0,1};
Point(4) = {L_AIR + L_Ag + L_AIR,0,0,1};



Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};


Physical Line("AIR1") = {1};
Physical Line("Ag") = {2};
Physical Line("AIR2") = {3};




Physical Point("point1") = {1};
Physical Point("point2") = {2};
Physical Point("point3") = {3};
Physical Point("point4") = {4};







