/*Point(1) = {-1,0,0,0.02};
Point(2) = {0,0,0,0.002};
Point(3) = {10,0,0,0.2};
Point(4) = {11,0,0,0.02};
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};


Physical Line("n_side") = {1};
Physical Line("p_side") = {2};
Physical Line("p_contactlayer") = {3};
Physical Point("anode") = {4};
Physical Point("cathode") = {1};*/


L = 10;
L_2 = L / 2;

nn = 100;
dx = L / nn;

Point(1) = {-L/2,0,0,dx};
Point(2) = {0,0,0,dx};
Point(3) = {L/2,0,0,dx};

Line(1) = {1,2};
Line(2) = {2,3};

Physical Line("n_side") = {1};
Physical Line("p_side") = {2};
Physical Point("anode") = {3};
Physical Point("cathode") = {1};

