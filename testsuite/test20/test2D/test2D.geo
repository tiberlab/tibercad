// Gmsh project created on Thu Feb 21 2012
lc1 = 0.25;
lc2 = 0.35;

c = 7;
br = 3;
we = 10;
l = 50;

Point(1) = {0, 0, 0, lc1};
Point(2) = {c, 0, 0, lc1};
Point(3) = {c + br, 0, 0, lc1};
Point(4) = {c + br + we, 0, 0, lc1};
Point(5) = {c + br + we + br, 0, 0, lc1};
Point(6) = {c + br + we + br + c, 0, 0, lc1};
Point(7) = {-l, 0, 0, lc2};
Point(8) = {c + br + we + br + c + l, 0, 0, lc2}; 


Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};
Line(6) = {7,1};
Line(7) = {6,8};


Physical Surface("dd_source_region") = {41};
Physical Surface("source_region") = {33};
Physical Surface("barrier1") = {29};
Physical Surface("well") = {25};
Physical Surface("barrier2") = {21};
Physical Surface("drain_region") = {17};
Physical Surface("dd_drain_region") = {37};
Physical Surface("oxide1") = {65, 69, 73, 77, 81};
Physical Surface("oxide2") = {45, 49, 53, 57, 61};

Physical Line("source") = {40};
Physical Line("int_source") = {32, 80, 60};
Physical Line("int_drain") = {15, 63, 43};
Physical Line("drain") = {35};

Extrude {0, 7, 0}{ 
 Line{1, 2, 3 ,4, 5, 6, 7}; 
}

Extrude {0, -3, 0}{
 Line{1, 2, 3, 4, 5};
}

Extrude {0, 3, 0}{
 Line{14, 18, 22, 26, 30};
}
