// the box base side length
Lbox = 20; 

// the box height
Hbox = 100; 

// the intrinsic layer thickness
Hintr = 20; 

// char. length at contacts
lc = 4;

// char. length at intrinsic layer boundary
li = 3;

c = Lbox / 2;
d = Hintr / 2;
e = Hbox / 2;


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
Surface Loop(117) = {108, 110, 112, 114, 106, 80};
Volume(118) = {117};
Surface Loop(121) = {98, 120, 100, 104, 102, 106};
Volume(122) = {121};


Physical Volume("nside") = {116};
Physical Volume("pside") = {122};
Physical Volume("intrinsic") = {118};

Physical Surface("cathode") = {70};
Physical Surface("anode") = {98};
Physical Surface("int_cathode") = {80};
Physical Surface("int_anode") = {106};


