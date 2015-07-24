lsub=0.03;
lacc=0.002;
lct=0.0005;
lg=0.0015;
lh=0.005;
lc=0.0005;

Lg_2 = 0.0375;
d = 0.01;
Ls = 0.1;
h = 0.25;
b = 0.0025;
o = 0.005;

xd = Lg_2 + d;
xd2 = Lg_2 + d / 2;
xmax = xd + Ls - d;

Point(1) = {0, -h, 0, lsub};
Point(2) = {0, 0, 0, lc};
Point(3) = {xmax,-h,0.0,lsub};
Point(4) = {-xmax,-h,0.0,lsub};
Point(5) = {xmax,0,0.0,lh};
Point(6) = {-xmax,0,0.0,lh};
Point(7) = {-xd,0,0.0,lct};
Point(8) = {-Lg_2 - o,0,0.0,lc};
Point(9) = {Lg_2 + o,0,0.0,lc};
Point(10) = {xd,0,0.0,lct};
Point(11) = {-Lg_2,0,0,lc};
Point(12) = {Lg_2,0,0,lc};
Point(13) = {xmax,-0.02,0,lh};
Point(14) = {-xmax,-0.02,0,lh};
Point(15) = {-xd2,-0.02,0.0,lacc};
Point(16) = {xd2,-0.02,0.0,lacc};
Point(17) = {-Lg_2 - o,b,0,lg};
Point(18) = {Lg_2 + o,b,0,lg};
Point(19) = {0, b, 0, lg};

Line(1) = {4,1};
Line(2) = {3,13};
Line(6) = {4,14};
Line(7) = {10,9};
Line(8) = {12,2};
Line(9) = {8,7};
Line(10) = {11,8};
Line(11) = {9,12};
Line(13) = {7,6};
Line(19) = {5,10};
Line(28) = {1,3};
Line(29) = {2,11};
Line(30) = {14,15};
Line(31) = {15,11};
Line(32) = {14,6};
Line(33) = {16,12};
Line(34) = {16,13};
Line(35) = {13,5};
Line(36) = {8,17};
Line(37) = {9,18};
Line(38) = {18,19};
Line(39) = {19,17};




Line Loop(40) = {28,2,-34,33,8,29,-31,-30,-6,1};
Plane Surface(41) = {40};
Line Loop(42) = {30,31,9,10,13,-32};
Plane Surface(47) = {42};
Line Loop(43) = {34,35,7,11,19,-33};
Plane Surface(44) = {43};
Line Loop(45) = {8,29,36,10,-39,-38,-37,11};
Plane Surface(46) = {45};

Physical Surface("substrate") = {41}; // n-Si
Physical Surface("contact") = {44,47}; // n+-Si
Physical Surface("oxide") = {46}; // SiO2
Physical Line("source") = {13}; // source
Physical Line("gate") = {39,38}; // gate
Physical Line("drain") = {19}; // drain
Physical Line("backcontact") = {1,28}; // drain
