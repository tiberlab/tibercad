a0 = 2;
a = 10;
a1 = 5;
a2 = 15;

d = 1;
z0 = 0;
h2 = 4;
h = 10;
h1 = 4;
h2 = 4;
h3 = 40;
N = 10;
N0 = 3;
N1 = 4;

//p7_0--------p7----------l9---------p8
//|            |                      |     
//|            |                      |
//|            |                      |
//p3_0-------- p3 ---p2               |
//|            |      \               |
//|            |       \              |
//|            |        \             |
//|            |         \            |
//|            |          \           |
//p0_0--------p0----------p1---------p4
//|            |                      |
//p5_0--------p5--------l6------------p6     
//
//
//
//

p0_0 = newp;
Point(p0_0) = {0, 0, z0, d};

p3_0 = newp;
Point(p3_0) = {0, 0, z0+h, d};

p5_0 = newp;
Point(p5_0) = {0, 0, z0-h1, d};

p7_0 = newp;
Point(p7_0) = {0, 0, z0+h+h2, d};


l1_0 = newl;
Line(l1_0) = {p0_0, p3_0};


l7_0 = newl;
Line(l7_0) = {p5_0, p0_0};


l8_0 = newl; 
Line(l8_0) = {p3_0, p7_0};

t1[] = Extrude {a0,0,0} {Line{l1_0}; Layers{N0}; Recombine; };

s1_0 = t1[1];

l1 = t1[0];



p0 = newp - 2;

p3 = newp - 1;

t7[] = Extrude {a0,0,0} {Line{l7_0}; Layers{N0}; Recombine; };

s2_0 = t7[1];
l7 = t7[0];
l7_1 = l7 +1;

p5 = newp - 1;

t8[] = Extrude {a0,0,0} {Line{l8_0}; Layers{N0}; Recombine; };

s3_0 = t8[1];
l8 = t8[0];

p7 = newp - 1;



p1 = newp;
Point(p1) = {a, 0, z0, d};


p2 = newp;
Point(p2) = {a1, 0, z0+h, d};


p4 = newp;
Point(p4) = {a2, 0, z0, d};


p6 = newp;
Point(p6) = {a2, 0,  z0-h1, d};


p8 = newp;
Point(p8) = {a2, 0, z0+h+h2, d};


l2 = newl;
Line(l2) = {p3, p2};
l3 = newl;
Line(l3) = {p2, p1};
l4 = newl;
Line(l4) = {p1, p4};

l5 = newl;
Line(l5) = {p4, p6};

l6 = newl;
Line(l6) = {p6, p5};


loop = newreg;


Line Loop(loop) = {l1, l2, l3, l4, l5, l6, l7};

s1 = news;
Plane Surface(s1) = {loop};




l9 = newl;
Line(l9) = {p7, p8};


l10 = newl;
Line(l10) = {p8, p4};



loop1 = newreg;
Line Loop(loop1) = {l8, l9, l10, -l4, -l3, -l2};
s2 = news;
Plane Surface(s2) = {loop1};

t1[] = Extrude {{0,0,1}, {0,0,0}, 2*Pi/3} {  Surface{s1};  Layers{N}; Recombine; };
t2[] = Extrude {{0,0,1}, {0,0,0}, 2*Pi/3} { Surface{t1[0]}; Layers{N}; Recombine; };
t3[] = Extrude {{0,0,1}, {0,0,0}, 2*Pi/3} { Surface{t2[0]}; Layers{N}; Recombine; };
vol_pyr1[] = {t1[1],t2[1], t3[1] };

t1[] = Extrude {{0,0,1}, {0,0,0}, 2*Pi/3} {  Surface{s1_0};  Layers{N}; Recombine; };
t2[] = Extrude {{0,0,1}, {0,0,0}, 2*Pi/3} { Surface{t1[0]}; Layers{N}; Recombine; };
t3[] = Extrude {{0,0,1}, {0,0,0}, 2*Pi/3} { Surface{t2[0]}; Layers{N}; Recombine; };
vol_pyr2[] = {t1[1],t2[1], t3[1] };

t1[] = Extrude {{0,0,1}, {0,0,0}, 2*Pi/3} {  Surface{s2_0};  Layers{N}; Recombine; };
t2[] = Extrude {{0,0,1}, {0,0,0}, 2*Pi/3} { Surface{t1[0]}; Layers{N}; Recombine; };
t3[] = Extrude {{0,0,1}, {0,0,0}, 2*Pi/3} { Surface{t2[0]}; Layers{N}; Recombine; };
vol_pyr3[] = {t1[1],t2[1], t3[1] };


volume_pyramid[] = {vol_pyr1[], vol_pyr2[], vol_pyr3[]}; 

Physical Volume(1) = {volume_pyramid[]};

t1[] = Extrude {{0,0,1}, {0,0,0}, 2*Pi/3} {  Surface{s3_0};  Layers{N}; Recombine; };
t2[] = Extrude {{0,0,1}, {0,0,0}, 2*Pi/3} { Surface{t1[0]}; Layers{N}; Recombine; };
t3[] = Extrude {{0,0,1}, {0,0,0}, 2*Pi/3} { Surface{t2[0]}; Layers{N}; Recombine; };
vol_c1[] = {t1[1],t2[1],t3[1]}; 

t1[] = Extrude {{0,0,1}, {0,0,0}, 2*Pi/3} {  Surface{s2};  Layers{N}; Recombine; };
t2[] = Extrude {{0,0,1}, {0,0,0}, 2*Pi/3} { Surface{t1[0]}; Layers{N}; Recombine; };
t3[] = Extrude {{0,0,1}, {0,0,0}, 2*Pi/3} { Surface{t2[0]}; Layers{N}; Recombine; };
vol_c2[] = {t1[1],t2[1],t3[1]};




t[] = Extrude {0,0,-h2} {Line{l6}; Layers{N1}; Recombine; };

s3 = t[1];

p[] = Extrude {0,0,-h3} {Line{t[0]}; Layers{N1}; Recombine; };

s4 = p[1];


t[] = Extrude {0,0,-h2} {Line{l7_1}; Layers{N1}; Recombine; };

s5 = t[1];

p[] = Extrude {0,0,-h3} {Line{t[0]}; Layers{N1}; Recombine; };

s6 = p[1];

t1[] = Extrude {{0,0,1}, {0,0,0}, 2*Pi/3} {  Surface{s3};  Layers{N}; Recombine; };
t2[] = Extrude {{0,0,1}, {0,0,0}, 2*Pi/3} { Surface{t1[0]}; Layers{N}; Recombine; };
t3[] = Extrude {{0,0,1}, {0,0,0}, 2*Pi/3} { Surface{t2[0]}; Layers{N}; Recombine; };
vol_c3[] = {t1[1],t2[1],t3[1]};


t1[] = Extrude {{0,0,1}, {0,0,0}, 2*Pi/3} {  Surface{s5};  Layers{N}; Recombine; };
t2[] = Extrude {{0,0,1}, {0,0,0}, 2*Pi/3} { Surface{t1[0]}; Layers{N}; Recombine; };
t3[] = Extrude {{0,0,1}, {0,0,0}, 2*Pi/3} { Surface{t2[0]}; Layers{N}; Recombine; };
vol_c4[] = {t1[1],t2[1],t3[1]};

Physical Volume(2) = {vol_c1[],vol_c2[],vol_c3[], vol_c4[]};


t1[] = Extrude {{0,0,1}, {0,0,0}, 2*Pi/3} {  Surface{s4};  Layers{N}; Recombine; };
t2[] = Extrude {{0,0,1}, {0,0,0}, 2*Pi/3} { Surface{t1[0]}; Layers{N}; Recombine; };
t3[] = Extrude {{0,0,1}, {0,0,0}, 2*Pi/3} { Surface{t2[0]}; Layers{N}; Recombine; };
vol_c3[] = {t1[1],t2[1],t3[1]};


t1[] = Extrude {{0,0,1}, {0,0,0}, 2*Pi/3} {  Surface{s6};  Layers{N}; Recombine; };
t2[] = Extrude {{0,0,1}, {0,0,0}, 2*Pi/3} { Surface{t1[0]}; Layers{N}; Recombine; };
t3[] = Extrude {{0,0,1}, {0,0,0}, 2*Pi/3} { Surface{t2[0]}; Layers{N}; Recombine; };
vol_c4[] = {t1[1],t2[1],t3[1]};

Physical Volume(3) = {vol_c3[], vol_c4[]};


p[] = Extrude {0,0,h3} {Line{l9}; Layers{N1}; Recombine; };

s4 = p[1];

t1[] = Extrude {{0,0,1}, {0,0,0}, 2*Pi/3} {  Surface{s4};  Layers{N}; Recombine; };
t2[] = Extrude {{0,0,1}, {0,0,0}, 2*Pi/3} { Surface{t1[0]}; Layers{N}; Recombine; };
t3[] = Extrude {{0,0,1}, {0,0,0}, 2*Pi/3} { Surface{t2[0]}; Layers{N}; Recombine; };
vol_c3[] = {t1[1],t2[1],t3[1]};


p[] = Extrude {0,0,h3} {Line{14}; Layers{N1}; Recombine; };

s4 = p[1];

t1[] = Extrude {{0,0,1}, {0,0,0}, 2*Pi/3} {  Surface{s4};  Layers{N}; Recombine; };
t2[] = Extrude {{0,0,1}, {0,0,0}, 2*Pi/3} { Surface{t1[0]}; Layers{N}; Recombine; };
t3[] = Extrude {{0,0,1}, {0,0,0}, 2*Pi/3} { Surface{t2[0]}; Layers{N}; Recombine; };
vol_c4[] = {t1[1],t2[1],t3[1]};

Physical Volume(4) = {vol_c3[], vol_c4[]};

Physical Surface(1) = {513, 535, 557, 610, 593, 576};