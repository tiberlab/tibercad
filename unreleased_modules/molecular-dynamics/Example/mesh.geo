S= 5;
L1 = 20;
L2 = 20;

Ns = 3;
Nl = 10;

p = newp; Point(newp) = {0, 0, 0,1};
t[] = Extrude{S,0,0}{ Point{p};Layers{Ns};Recombine;};

t[] = Extrude{0,S,0}{ Line{t[1]};Layers{Ns};Recombine;};
Physical Surface("Left") = {t[1]};

t[] = Extrude{0,0,L1}{ Surface{t[1]};Layers{Nl};Recombine;};
Physical Volume("Bulk1") = {t[1]};

t[] = Extrude{0,0,L2}{ Surface{t[0]};Layers{Nl};Recombine;};
Physical Volume("Bulk2") = {t[1]};

t[] = Extrude{0,0,L2}{ Surface{t[0]};Layers{Nl};Recombine;};
Physical Volume("Bulk3") = {t[1]};

t[] = Extrude{0,0,L1}{ Surface{t[0]};Layers{Nl};Recombine;};
Physical Volume("Bulk4") = {t[1]};

Physical Surface("Right") = {t[0]};
