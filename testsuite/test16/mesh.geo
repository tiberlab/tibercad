R = 0.150 * 1;
H = 4 * 1;

M_y = 2 * H;
M_x = 10 * R;
A_y = 2 * M_y; //2
A_x = 10 * M_x; //10



N = 0.05;
N2 = 0.01;
M = 0.1; 
A = 3;

p0 = newp; Point(p0) = {0,0,0,N};
p1 = newp; Point(p1) = {R,0,0,N};
p2 = newp; Point(p2) = {R,H,0,N2};
p3 = newp; Point(p3) = {0,H,0,N2};

p4 = newp; Point(p4) = {M_x,0,0,M};
p5 = newp; Point(p5) = {M_x,M_y,0,M};
p6 = newp; Point(p6) = {0,M_y,0,M};

p7 = newp; Point(p7) = {A_x,0,0,A};
p8 = newp; Point(p8) = {A_x,A_y,0,A};
p9 = newp; Point(p9) = {0,A_y,0,A};


Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(5) = {2, 5};
Line(6) = {5, 8};
Line(7) = {8, 9};
Line(8) = {9, 10};
Line(9) = {10, 7};
Line(10) = {7, 4};
Line(11) = {7, 6};
Line(12) = {6, 5};
Line Loop(13) = {11, 12, -5, 2, 3, -10};
Plane Surface(14) = {13};
Line Loop(15) = {4, 1, 2, 3};
Plane Surface(16) = {15};
Line Loop(17) = {9, 11, 12, 6, 7, 8};
Plane Surface(18) = {17};
Physical Line("Lower") = {1};
Physical Line("Upper") = {3};
Physical Surface("Column") = {16};
#Physical Surface("Air") = {14, 18};
