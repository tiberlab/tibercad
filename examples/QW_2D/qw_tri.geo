ymax = 51.5;
width = 20.0;
well = 3.0;
lc = width / 5.0;
lw = width / 40.0;

xmax = width / 2;
xmin = -xmax;


Point(1) = {xmin, -ymax, 0, lc};
Point(2) = {xmax, -ymax, 0, lc};
Point(3) = {xmin, -well / 2, 0, lw};
Point(4) = {xmax, -well / 2, 0, lw};
Point(5) = {xmin, well / 2, 0, lw};
Point(6) = {xmax, well / 2, 0, lw};
Point(7) = {xmax, ymax, 0, lc};
Point(8) = {xmin, ymax, 0, lc};

Line(5) = {1,2};
Line(6) = {2,4};
Line(7) = {3,4};
Line(8) = {3,1};
Line(9) = {4,6};
Line(10) = {5,6};
Line(11) = {5,3};
Line(12) = {6,7};
Line(13) = {7,8};
Line(14) = {8,5};

Line Loop(15) = {5,6,-7,8};
Plane Surface(16) = {15};
Line Loop(17) = {7,9,-10,11};
Plane Surface(18) = {17};
Line Loop(19) = {10,12,13,14};
Plane Surface(20) = {19};

Physical Surface(1) = {16};
Physical Surface(2) = {18};
Physical Surface(3) = {20};
Physical Line(1) = {13};
Physical Line(2) = {5};
