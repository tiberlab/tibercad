// the box base side length
Lbox = 10; // 50;

// the box height
Hbox = 20; //100;

// the intrinsic layer thickness
Hintr = 6; //30;

// the quantum cube side length
Lqbox = 4.0; // 5;  //15;

// the radius of the ball
R = 1.0;  // 1.5;

// char. length at contacts
lc = 5;

// char. length for quantum box
lq = 1;

// char. length for quantum ball
lb = 0.15;

// char. length at intrinsic layer boundary
li = 3;

b = Lqbox / 2;
c = Lbox / 2;
d = Hintr / 2;
e = Hbox / 2;

Point (1) = {0, 0, 0, lb};
Point (2) = {R, 0, 0, lb};
Point (3) = {0, R, 0, lb};
Point (4) = {-R, 0, 0, lb};
Point (5) = {0, -R, 0, lb};
Point (6) = {0, 0, -R, lb};
Point (7) = {0, 0, R, lb};
Point (9) = {-b, -b, -b, lq};
Point (10) = {b, -b, -b, lq};
Point (11) = {b, b, -b, lq};
Point (12) = {-b, b, -b, lq};
Point (13) = {-b, b, b, lq};
Point (14) = {-b, -b, b, lq};
Point (15) = {b, -b, b, lq};
Point (16) = {b, b, b, lq};
Circle (1) = {2, 1, 3};
Circle (2) = {3, 1, 4};
Circle (3) = {4, 1, 5};
Circle (4) = {5, 1, 2};
Circle (5) = {3, 1, 6};
Circle (6) = {6, 1, 5};
Circle (7) = {5, 1, 7};
Circle (8) = {7, 1, 3};
Circle (9) = {2, 1, 7};
Circle (10) = {7, 1, 4};
Circle (11) = {4, 1, 6};
Circle (12) = {6, 1, 2};
Line (31) = {14, 15};
Line (32) = {15, 16};
Line (33) = {16, 13};
Line (34) = {13, 14};
Line (35) = {10, 11};
Line (36) = {11, 12};
Line (37) = {12, 9};
Line (38) = {9, 10};
Line (39) = {10, 15};
Line (40) = {16, 11};
Line (41) = {12, 13};
Line (42) = {14, 9};
Line Loop (14) = {2, -10, 8};
Ruled Surface (14) = {14};
Line Loop (16) = {10, 3, 7};
Ruled Surface (16) = {16};
Line Loop (18) = {-8, -9, 1};
Ruled Surface (18) = {18};
Line Loop (20) = {-11, -2, 5};
Ruled Surface (20) = {20};
Line Loop (22) = {-5, -1, -12};
Ruled Surface (22) = {22};
Line Loop (24) = {-3, 11, 6};
Ruled Surface (24) = {24};
Line Loop (26) = {-7, 4, 9};
Ruled Surface (26) = {26};
Line Loop (28) = {-4, -6, 12};
Ruled Surface (28) = {28};
Line Loop (44) = {39, 32, 40, -35};
Plane Surface (44) = {44};
Line Loop (46) = {38, 35, 36, 37};
Plane Surface (46) = {46};
Line Loop (48) = {39, -31, 42, 38};
Plane Surface (48) = {48};
Line Loop (50) = {34, 42, -37, 41};
Plane Surface (50) = {50};
Line Loop (52) = {33, -41, -36, -40};
Plane Surface (52) = {52};
Surface Loop (30) = {28, 26, 16, 14, 20, 24, 22, 18};
Volume (30) = {30};
Line Loop(53) = {32,33,34,31};
Plane Surface(54) = {53};
Surface Loop(55) = {44,48,54,52,50,46};
Volume(56) = {55,30};


Point(17) = {-c, -c, -e, lc};
Point(18) = {c, -c, -e, lc};
Point(19) = {-c, c, -e, lc};
Point(20) = {c, c, -e, lc};
Point(21) = {-c, -c, -d, li};
Point(22) = {c, -c, -d, li};
Point(23) = {-c, c, -d, li};
Point(24) = {c, c, -d, li};
Line(57) = {17, 18};
Line(58) = {18, 20};
Line(59) = {20, 19};
Line(60) = {19, 17};
Line(61) = {17, 21};
Line(62) = {21, 23};
Line(63) = {23, 19};
Line(64) = {21, 22};
Line(65) = {22, 18};
Line(66) = {22, 24};
Line(67) = {20, 24};
Line(68) = {24, 23};
Line Loop(69) = {57, 58, 59, 60};
Plane Surface(70) = {69};
Line Loop(71) = {61, 64, 65, -57};
Plane Surface(72) = {71};
Line Loop(73) = {67, -66, 65, 58};
Plane Surface(74) = {73};
Line Loop(75) = {61, 62, 63, 60};
Plane Surface(76) = {75};
Line Loop(77) = {63, -59, 67, 68};
Plane Surface(78) = {77};
Line Loop(79) = {64, 66, 68, -62};
Plane Surface(80) = {79};

Point(31) = {-c, -c, e, lc};
Point(32) = {c, -c, e, lc};
Point(33) = {-c, c, e, lc};
Point(34) = {c, c, e, lc};
Point(35) = {-c, -c, d, li};
Point(36) = {c, -c, d, li};
Point(37) = {-c, c, d, li};
Point(38) = {c, c, d, li};
Line(81) = {32, 34};
Line(82) = {34, 33};
Line(83) = {33, 31};
Line(84) = {31, 32};
Line(85) = {32, 36};
Line(86) = {34, 38};
Line(87) = {38, 36};
Line(88) = {31, 35};
Line(89) = {35, 37};
Line(90) = {37, 33};
Line(91) = {37, 38};
Line(92) = {36, 35};
Line(93) = {35, 21};
Line(94) = {23, 37};
Line(95) = {38, 24};
Line(96) = {36, 22};
Line Loop(97) = {82, 83, 84, 81};
Plane Surface(98) = {97};
Line Loop(99) = {81, 86, 87, -85};
Plane Surface(100) = {99};
Line Loop(101) = {83, 88, 89, 90};
Plane Surface(102) = {101};
Line Loop(103) = {88, -92, -85, -84};
Plane Surface(104) = {103};
Line Loop(105) = {91, 87, 92, 89};
Plane Surface(106) = {105};
Line Loop(107) = {96, -64, -93, -92};
Plane Surface(108) = {107};
Line Loop(109) = {87, 96, 66, -95};
Plane Surface(110) = {109};
Line Loop(111) = {95, 68, 94, 91};
Plane Surface(112) = {111};
Line Loop(113) = {94, -89, 93, 62};
Plane Surface(114) = {113};
Line Loop(119) = {86, -91, 90, -82};
Plane Surface(120) = {119};

Surface Loop(115) = {72, 76, 78, 70, 74, 80};
Volume(116) = {115};
Surface Loop(117) = {108, 110, 112, 114, 106, 48, 44, 54, 52, 50, 46, 80};
Volume(118) = {117};
Surface Loop(121) = {98, 120, 100, 104, 102, 106};
Volume(122) = {121};


Physical Volume("ball") = {30};
Physical Volume("qbox") = {56};
Physical Volume("nside") = {116};
Physical Volume("pside") = {122};
Physical Volume("intrinsic") = {118};

Physical Surface("cathode") = {70};
Physical Surface("anode") = {98};
