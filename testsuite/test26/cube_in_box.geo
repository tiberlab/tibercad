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

lx = (na-3)*a/2;
ly = (nb-3)*b/2;
lz = (nc-3)*c/2;


// char. length for box
lq = 0.50;

// char. length for cube
lb = 0.50;

Point (19) = {-lx, -ly, -lz, lb};
Point (20) = {lx, -ly, -lz, lb};
Point (21) = {lx, ly, -lz, lb};
Point (22) = {-lx, ly, -lz, lb};
Point (23) = {-lx, ly , lz, lb};
Point (24) = {-lx, -ly, lz, lb};
Point (25) = {lx, -ly, lz, lb};
Point (26) = {lx, ly, lz, lb};

Line (61) = {24, 25};
Line (62) = {25, 26};
Line (63) = {26, 23};
Line (64) = {23, 24};
Line (65) = {20, 21};
Line (66) = {21, 22};
Line (67) = {22, 19};
Line (68) = {19, 20};
Line (69) = {20, 25};
Line (70) = {26, 21};
Line (71) = {22, 23};
Line (72) = {24, 19};
// 6 faces of the small cube
Line Loop (84) = {69, 62, 70, -65};
Plane Surface (84) = {84};
Line Loop (86) = {68, 65, 66, 67};
Plane Surface (86) = {86};
Line Loop (88) = {69, -61, 72, 68};
Plane Surface (88) = {88};
Line Loop (90) = {64, 72, -67, 71};
Plane Surface (90) = {90};
Line Loop (92) = {63, -71, -66, -70};
Plane Surface (92) = {92};
Line Loop (93) = {62,63,64,61};
Plane Surface (94) = {93};

Surface Loop(130) = {84,88,94,92,90,86};
Volume (130) = {130};

// BOX definition
Point (9) = {-Lx, -Ly, -Lz, lq};
Point (10) = {Lx, -Ly, -Lz, lq};
Point (11) = {Lx, Ly, -Lz, lq};
Point (12) = {-Lx, Ly, -Lz, lq};
Point (13) = {-Lx, Ly , Lz, lq};
Point (14) = {-Lx, -Ly, Lz, lq};
Point (15) = {Lx, -Ly, Lz, lq};
Point (16) = {Lx, Ly, Lz, lq};
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

// 6 faces of the box
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
Volume(56) = {55,130};

Physical Volume("cube") = {130};
Physical Volume("box") = {56};


