L_AIR = 20;
L_Glass = 2e6;
dummy_point = L_Glass/2;



Point(1) = {0,0,0,1};
Point(2) = {L_AIR,0,0,1};
Point(3) = {dummy_point + L_AIR,0,0,1e4};
Point(4) = {L_Glass + L_AIR,0,0,1};
Point(5) = {L_AIR + L_Glass + L_AIR,0,0,1};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};

Physical Line("AIR") = {1,4};
Physical Line("Glass") = {2,3};

Physical Point("point1") = {1};









