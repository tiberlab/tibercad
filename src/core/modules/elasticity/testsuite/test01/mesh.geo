
d = 50;
HC = 100;
HQ = 10;

//---------------------
Np = 7;
NC =40;
NQ = 20;
//-----------------

//1D
p = newp; Point(p)  = (-d/2,0,0);
t[] = Extrude {d/2,0,0} {Point{p};  Layers{Np};Recombine;};
t[] = Extrude {0,d/2,0} {Line{t[1]};  Layers{Np};Recombine;};
Physical Surface("Base")={t[1]};

t[] = Extrude {0,0,HC} {Surface{t[1]};  Layers{NC};Recombine;};
V[0] = t[1];
t[] = Extrude {0,0,HQ} {Surface{t[0]};  Layers{NQ};Recombine;};
Q = t[1];
t[] = Extrude {0,0,HC} {Surface{t[0]};  Layers{NC};Recombine;};
V[1] = t[1];

Physical Volume("Contact")={V[]};
Physical Volume("Well")={Q};







