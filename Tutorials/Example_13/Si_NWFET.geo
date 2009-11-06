// Gmsh project created on Fri Jun 26 09:47:10 2009
lc1 =3;
lc2 = 3;
lc3 = 15;
Point(1) = {10, 10, 0, lc1};
Point(2) = {0, 0, 0, lc1};
Point(3) = {20, 0, 0, lc1};
Point(4) = {10, -10, 0, lc1};
Point(5) = {10, -11.6, 0, lc1};
Point(6) = {10, 11.6, 0, lc2};
Point(7) = {-1.6, 0, 0, lc2};
Point(8) = {21.6, 0, 0, lc2};
Point(9) = {10, 0, 0, lc3};
Circle(1) = {5, 9, 8};
Circle(2) = {8, 9, 6};
Circle(3) = {6, 9, 7};
Circle(4) = {7, 9, 5};
Circle(5) = {4, 9, 3};
Circle(6) = {3, 9, 1};
Circle(7) = {1, 9, 2};
Circle(8) = {2, 9, 4};
Extrude {0, 0, 60} {
  Point{3, 8, 9};
}
Extrude {0, 0, 60} {
  Point{1, 6, 7, 2};
}
Extrude {0, 0, 60} {
  Point{5, 4};
}
Extrude {0, 0, 60} {
  Line{2, 6, 7, 3, 4, 8, 5, 1};
}
Line(50) = {6, 1};
Line(51) = {4, 5};
Extrude {0, 0, 60} {
  Line{50, 51};
}
Line Loop(60) = {2, 50, -6, -5, 51, 1};
Plane Surface(61) = {60};
Line Loop(62) = {3, 4, -51, -8, -7, -50};
Plane Surface(63) = {62};
Line Loop(64) = {8, 5, 6, 7};
Plane Surface(65) = {64};
Extrude {0, 0, 60} {
  Surface{65};
}
Extrude {0, 0, 60} {
  Surface{61, 63};
}
Extrude {0, 0, -20} {
  Point{1, 3, 4, 2};
}
Extrude {0, 0, -20} {
  Line{6, 5, 8, 7};
}
Extrude {0, 0, -20} {
  Surface{65};
}
Extrude {0, 0, -20} {
  Point{9};
}
Extrude {0, 0, -30} {
  Point{24, 21, 20, 19, 22};
}
Extrude {0, 0, -30} {
  Line{156, 160, 164, 168};
}
Extrude {0, 0, -30} {
  Surface{193};
}
Extrude {0, 0, 20} {
  Point{16, 18, 10, 13, 12};
}
Extrude {0, 0, 20} {
  Line{38, 26, 22, 42};
}
Extrude {0, 0, 20} {
  Surface{87};
}
Extrude {0, 0, 30} {
  Point{31, 32, 33, 34, 30};
}
Extrude {0, 0, 30} {
  Line{247, 243, 255, 251};
}
Extrude {0, 0, 30} {
  Surface{280};
}
Physical Surface("drain") = {323};
Physical Surface("gate") = {49, 21, 33, 37};
Physical Surface("source") = {237};

Physical Volume("oxide") = {3, 2};
Physical Volume("channel") = {1, 4, 6};
Physical Volume("drain_region") = {7};
Physical Volume("source_region") = {5};
