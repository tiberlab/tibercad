lcc = 0.2;
lci = 0.2;
lcb = 0.1;
lbe = 0.1;
lec = 0.5;
leci = 0.2;
la = 1;
lbc = 0.2;

L = 17;
H = 7.5;
Hcd = 2;
Hcu = 4;
W = 0.5;
Lec = 4;
Le = 6;
Lb = 4;
He = H - Hcd - Hcu - W;

Point(1) = {0, 0, 0, lcc};
Point(2) = {L, 0, 0, lcc};
Point(3) = {0, Hcd, 0, lci};
Point(4) = {L, Hcd, 0, lci};
Point(5) = {L, H, 0, la};
Point(6) = {0, Hcd+Hcu, 0, lcb};
Point(7) = {0, Hcd+Hcu+W, 0, lbe};
Point(8) = {0, H, 0, lec};
Point(9) = {Lec, H, 0, leci};
Point(10) = {Le, H, 0, lbe};
Point(11) = {Le+Lb, H, 0, lcb};
Point(12) = {Le, H-(He-1), 0, lbe};
Point(13) = {Le+Lb, H-(He+W-1), 0, lcb};
Point(14) = {Le+1, H, 0, lbc};
Point(15) = {Le+Lb-1, H, 0, lbc};

Line(1) = {1, 2};
Line(2) = {2, 4};
Line(3) = {4, 3};
Line(4) = {3, 1};
Line(5) = {4, 5};
Line(8) = {10, 9};
Line(9) = {9, 8};
Line(10) = {8, 7};
Line(11) = {7, 6};
Line(12) = {6, 3};
Line(13) = {10, 12};
Line(14) = {11, 13};

Extrude { {0, 0, 1}, {Le-1,H-(He-1), 0}, -Pi/2} {Point {12}; }
Extrude { {0, 0, 1}, {Le+Lb-1, H-(He+W-1), 0}, -Pi/2} { Point {13}; }

p = newp; Point(p) = {L-5, H, 0, lcc};

Line(17) = {10, 14};
Line(18) = {14, 15};
Line(19) = {15, 11};
Line(20) = {18, 6};
Line(21) = {16, 7};
l = newl; Line(l) = {5, p};
l = newl; Line(l) = {p, 11};

Line Loop(22) = {1, 2, 3, 4};
Plane Surface(23) = {22};
Line Loop(24) = {3, -12, -20, -16, -14, -l, -l+1, -5};
Plane Surface(25) = {24};
Line Loop(26) = {20, -11, -21, -15, -13, 17, 18, 19, 14, 16};
Plane Surface(27) = {26};
Line Loop(28) = {8, 9, 10, -21, -15, -13};
Plane Surface(29) = {28};

Physical Surface("coll_doped") = {23};
Physical Surface("coll") = {25};
Physical Surface("b_reg") = {27};
Physical Surface("e_reg") = {29};

Physical Line("emitter") = {9};
Physical Line("base") = {18};
Physical Line("collector") = {22};
