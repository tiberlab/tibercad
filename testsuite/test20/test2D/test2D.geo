// Gmsh project created on Thu Feb 21 2012
lc1 =0.25;

c = 10;
br = 2;
we = 5;


Point(1) = {0, 0, 0, lc1};
Point(2) = {c, 0, 0, lc1};
Point(3) = {c + br, 0, 0, lc1};
Point(4) = {c + br + we, 0, 0, lc1};
Point(5) = {c + br + we + br, 0, 0, lc1};
Point(6) = {c + br + we + br + c, 0, 0, lc1};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};

Physical Surface("source_region") = {13};
Physical Surface("barrier1") = {17};
Physical Surface("well") = {21};
Physical Surface("barrier2") = {25};
Physical Surface("drain_region") = {29};

Physical Line("source") = {11};
Physical Line("int_source") = {12};
Physical Line("int_drain") = {24};
Physical Line("drain") = {28};

Extrude {0, 5, 0}{ 
  Line{1, 2, 3 ,4, 5}; 
}
