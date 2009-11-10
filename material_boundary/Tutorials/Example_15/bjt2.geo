lcc = 0.5;
lci = 0.2;
lcb = 0.01;
lcbo = 0.1;
lbe = 0.005;
lbm = 0.05;
lbeo = 0.1;
lec = 0.5;
leci = 0.2;
la = 1;
lbc = 0.2;

L = 17;
He = 2;
Hcd = 2;
Hcu = 4;
W = 0.5;
Lec = 4;
Le = 6;
Lb = 4;
H = He + Hcd + Hcu + W;
dc = 0.75;

Point(1) = {0, 0, 0, lcc};
Point(2) = {L, 0, 0, lcc};
Point(3) = {0, Hcd, 0, lci};
Point(4) = {L, Hcd, 0, lci};
Point(5) = {L, H, 0, la};
Point(6) = {0, Hcd+Hcu, 0, lcb};
Point(7) = {0, Hcd+Hcu+W, 0, lbe};
Point(8) = {0, H, 0, lec};
Point(9) = {Lec, H, 0, leci};
Point(11) = {Le+Lb, H, 0, lec};
Point(12) = {Le, H-(He-1), 0, lbe};
Point(14) = {Le+1, H, 0, lbc};
Point(16) = {0, Hcd+Hcu-dc, 0, lcbo};
Point(17) = {0, Hcd+Hcu+W+dc, 0, lbeo};
Point(18) = {0, Hcd+Hcu+W/2, 0, lbm};
Point(19) = {Le+Lb, Hcd, 0, lci};

Line(7) = {19, 3};
Line(8) = {1, 2};
Line(9) = {2, 4};
//Line(10) = {4, 19};
Line(11) = {3, 1};
Line(12) = {3, 16};
Line(13) = {16, 6};
Line(14) = {6, 18};
Line(15) = {18, 7};
Line(16) = {7, 17};
Line(17) = {17, 8};
Line(18) = {8, 9};




Extrude {Le-1, 0, 0} { Line {13,14,15,16}; Layers {20}; }
Extrude { {0, 0, 1}, {Le-1,H-(He-1), 0}, Pi/2} {Line {27,31}; Layers {10}; }
Extrude {0, He-1, 0} { Line {35,39}; Layers {10}; }
Extrude {Lb+1, 0, 0} { Line {19,23}; Layers {20}; }
Extrude { {0, 0, 1}, {Le+Lb+Sin(Pi/3),Hcd+Cos(Pi/3), 0}, 2*Pi/3} {Point {19}; }
Extrude {0, H-Hcd, 0} { Point {44}; }

p = newp; Point(p) = {L-5, H, 0, lcc};
Line(61) = {4, 5};



Line(62) = {5, 47};
Line(63) = {47, 46};
Line(64) = {19, 40};
Line(65) = {43, 11};
Line(66) = {11, 14};
Line(67) = {14, 36};
Line(68) = {39, 9};
Line(69) = {46, 11};
Line(70) = {44, 4};
Line Loop(71) = {8, 9, -70, -59, 7, 11};
Plane Surface(72) = {71};

Line Loop(73) = {70, 61, 62, 63, -60};
Plane Surface(74) = {73};
Line Loop(75) = {64, -52, -20, -12, -7};
Plane Surface(76) = {75};

Line Loop(77) = {65, 66, 67, -44, -36, 57};
Plane Surface(78) = {77};

Line Loop(79) = {49, 68, -18, -17, 33, 41};
Plane Surface(80) = {79};
Line Loop(81) = {60, 69, -65, -55, -51, -64, 59};
Plane Surface(82) = {81};



Physical Surface("coll_doped") = {72};
Physical Surface("coll") = {74,76,54,22};
Physical Surface("b_reg") = {78,58,46,38,30,26};
Physical Surface("e_reg") = {80,50,42,34};
Physical Surface("oxide") = {82};

Physical Line("emitter") = {18};
Physical Line("base") = {66};
Physical Line("collector") = {62};
