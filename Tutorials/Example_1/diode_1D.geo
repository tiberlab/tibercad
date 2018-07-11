/*Point(1) = {-25,0,0,0.1};
Point(2) = {0,0,0,0.05};
Point(3) = {25,0,0,0.1};
Line(1) = {1,2};
Line(2) = {2,3};


Physical Line("p_side") = {1};
Physical Line("n_side") = {2};
Physical Point("anode") = {1};
Physical Point("cathode") = {3};*/

/*L = 10;
d = 0.05;
dc = 0.01;


Point(1) = {-L/2, 0, 0, dx};
Point(2) = {0, 0, 0, dx};
Point(3) = {L/2, 0, 0, dx};

Line(1) = {1,2};
Line(2) = {2,3};

//t[] = Extrude{0,L,0}{Line{1,2}; Layers{d};};

Physical Surface("p_side") = {6};
Physical Surface("n_side") = {10};
Physical Line("anode") = {4};
Physical Line("cathode") = {9};*/

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
