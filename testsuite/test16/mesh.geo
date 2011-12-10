R = 0.150;
H = 4.0;
m = 0.01;
m1 = 0.01;
B = 3.7;
N = 10;
N2 = 20;


p0 = newp; Point(p0) = {0,0,0,m1};
t0[] = Extrude {R,0,0} { Point{p0}; Layers{N}; Recombine; };
t1[] = Extrude {0,0,B} { Line{t0[1]}; Layers{N2}; Recombine; };
t2[] = Extrude {0,0,H-B} { Line{t1[0]}; Layers{N2}; Recombine; };


t1[] = Extrude {{0,0,1}, {0,0,0}, 2*Pi/3} { Surface{t1[1],t2[1]}; Layers{N}; Recombine; };
t2[] = Extrude {{0,0,1}, {0,0,0}, 2*Pi/3} { Surface{43,26}; Layers{N}; Recombine; };
t2[] = Extrude {{0,0,1}, {0,0,0}, 2*Pi/3} { Surface{60,77}; Layers{N}; Recombine; };
Physical Surface("Lower") = {17,68,100};
Physical Surface("Upper") = {41,58,92};
Physical Volume("Column") = {1,2,3,4,5,6};
