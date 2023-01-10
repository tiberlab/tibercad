L_AIR = 100;
L_ITO = 200;
L_Si = 1e5;



lc = 5;
lc2 = 20;

x = -(L_AIR + L_ITO);
Point(1) = {x,0,0,lc};
Point(2) = {x + L_AIR,0,0,lc};
Point(3) = {x + L_ITO + L_AIR,0,0,lc};
Point(4) = {x + L_Si + L_ITO + L_AIR,0,0,lc2};
Point(5) = {x + L_AIR + L_Si + L_ITO + L_AIR,0,0,lc2};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};

Physical Line("AIR") = {1,4};
Physical Line("ITO") = {2};
Physical Line("Si") = {3};


Physical Point("front") = {1};
Physical Point("back") = {5};




