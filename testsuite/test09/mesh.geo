a = 20;

//0D
p = newp; Point(p[0]) = {-a/2  ,-a/2  ,0 ,0.1};
//1D
t[] = Extrude {0,a,0} {Point{p};  Layers{10};Recombine;};
//2D
t[] = Extrude {a,0,0} {Line{t[1]};  Layers{10};Recombine;};
Base = t[1];
//3D
t[] = Extrude {0,0,10} {Surface{t[1]};  Layers{10};Recombine;};
V[0]=t[1];
t[] = Extrude {0,0,10} {Surface{t[0]};  Layers{10};Recombine;};
Well=t[1];
t[] = Extrude {0,0,10} {Surface{t[0]};  Layers{10};Recombine;};
V[1]=t[1];

Physical Surface("Base")={Base};
Physical Volume("Bulk")={V[]};
Physical Volume("Well")={Well};






