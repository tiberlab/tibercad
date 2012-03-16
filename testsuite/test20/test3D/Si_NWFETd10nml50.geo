// Gmsh project created on Fri Jun 26 09:47:10 2009
lc1 =0.7;
lc2 = 1.7;
lc3 =0.35;
Point(1) = {5, 5, 0, lc1};
Point(2) = {0, 0, 0, lc1};
Point(3) = {10, 0, 0, lc1};
Point(4) = {5, -5, 0, lc1};
Point(5) = {5, -6, 0, lc2};
Point(6) = {5, 6, 0, lc2};
Point(7) = {-1, 0, 0, lc2};
Point(8) = {11, 0, 0, lc2};
Point(9) = {5, 0, 0, lc3};
Circle(1) = {5, 9, 8};
Circle(2) = {8, 9, 6};
Circle(3) = {6, 9, 7};
Circle(4) = {7, 9, 5};
Circle(5) = {4, 9, 3};
Circle(6) = {3, 9, 1};
Circle(7) = {1, 9, 2};
Circle(8) = {2, 9, 4};

Extrude {0, 0, 10} {
  Point{3, 8, 9};
}
Extrude {0, 0, 10} {
  Point{1, 6, 7, 2};
}
Extrude {0, 0, 10} {
  Point{5, 4};
}
Extrude {0, 0, 10} {
  Line{2, 6, 7, 3, 4, 8, 5, 1};
}

Line(50) = {6, 1};
Line(51) = {4, 5};

Extrude {0, 0, 8} {
  Line{50, 51};
}

Line Loop(60) = {2, 50, -6, -5, 51, 1};
Plane Surface(61) = {60};
Line Loop(62) = {3, 4, -51, -8, -7, -50};
Plane Surface(63) = {62};
Line Loop(64) = {8, 5, 6, 7};
Plane Surface(65) = {64};

Extrude {0, 0, 10} {
  Surface{65};
}
Extrude {0, 0, 10} {
  Surface{61, 63};
}
Extrude {0, 0, -4} {
  Point{1, 3, 4, 2};
}
Extrude {0, 0, -4} {
  Line{6, 5, 8, 7};
}
Extrude {0, 0, -4} {
  Surface{65};
}
Extrude {0, 0, -4} {
  Point{9};
}
Extrude {0, 0, -4} {
  Point{24, 21, 20, 19, 22};
}
Extrude {0, 0, -4} {
  Line{156, 160, 164, 168};
}
Extrude {0, 0, -4} {
  Surface{193};
}
Extrude {0, 0, 4} {
  Point{16, 18, 10, 13, 12};
}
Extrude {0, 0, 4} {
  Line{38, 26, 22, 42};
}
Extrude {0, 0, 4} {
  Surface{87};
}
Extrude {0, 0, 4} {
  Point{31, 32, 33, 34, 30};
}
Extrude {0, 0, 4} {
  Line{247, 243, 255, 251};
}
Extrude {0, 0, 4} {
  Surface{280};
}

Extrude {0, 0, 50} {
  Surface{318};
}

Extrude {0, 0, -50} {
  Surface{237};
}

Physical Surface("drain") = {318};
Physical Surface("gate") = {49, 21, 33, 37};
Physical Surface("source") = {237};
Physical Surface("inside_d") = {280};
Physical Surface("inside_s") = {193};
Physical Surface("ext_drain") = {340};
Physical Surface("ext_source") = {362};

Physical Volume("oxide") = {3, 2};
Physical Volume("channel1") = {1};
Physical Volume("channel4") = {4};
Physical Volume("channel6") = {6};
Physical Volume("drain_region") = {7};
Physical Volume("source_region") = {5};
Physical Volume("ext_source_region") = {9};
Physical Volume("ext_drain_region") = {8};
