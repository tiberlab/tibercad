//lattice constants of orthorhombic
a = 0.88556; 
b = 0.88556;
c = 1.26600;

// Ncells
na=8;
nb=8;
nc=6;

Lx = na*a/2;
Ly = nb*b/2;
Lz = nc*c/2;

// the radius of the ball
R = 2*Lx/3; 

// char. length for quantum box
lq = 0.50;

// char. length for quantum ball
lb = 0.50;

Point (1) = {0, 0, 0, lb};
Point (2) = {R, 0, 0, lb};
Point (3) = {0, R, 0, lb};
Point (4) = {-R, 0, 0, lb};
Point (5) = {0, -R, 0, lb};
Point (6) = {0, 0, -R, lb};
Point (7) = {0, 0, R, lb};
Point (9) = {-Lx, -Ly, -Lz, lq};
Point (10) = {Lx, -Ly, -Lz, lq};
Point (11) = {Lx, Ly, -Lz, lq};
Point (12) = {-Lx, Ly, -Lz, lq};
Point (13) = {-Lx, Ly , Lz, lq};
Point (14) = {-Lx, -Ly, Lz, lq};
Point (15) = {Lx, -Ly, Lz, lq};
Point (16) = {Lx, Ly, Lz, lq};
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

// 8 pieces of the sphere
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

Surface Loop (30) = {28, 26, 16, 14, 20, 24, 22, 18};
Volume (30) = {30};

// 6 faces of the cube
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
Line Loop (53) = {32,33,34,31};
Plane Surface (54) = {53};

Surface Loop(55) = {44,48,54,52,50,46};
//Volume(55) = {55};

Volume(56) = {55,30};


Physical Volume("ball") = {30};
Physical Volume("box") = {56};

