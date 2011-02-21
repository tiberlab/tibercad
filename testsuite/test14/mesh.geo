device_units = 1e-9;
mesh_units = 1e-9;
//---------
side    = 0.5;
lenght  = 5;

//---------------------
Nx = 10;
Ny = 10;
Nz = 20;
//-----------------


scale = device_units/mesh_units; 
a = side*scale;
len = lenght*scale; 


//0D
p = newp; Point(p[0]) = {-a/2  ,-a/2  , -len/2 ,0.1};
//1D
t[] = Extrude {0,a,0} {Point{p};  Layers{Nx};Recombine;};
prova = t[1];
//2D
t[] = Extrude {a,0,0} {Line{t[1]};  Layers{Ny};Recombine;};
s2 = t[1];
//3D
t[] = Extrude {0,0,len} {Surface{t[1]};  Layers{Nz};Recombine;};
V[0]=t[1];
s1 = t[0];
fs = t[5];

Physical Surface("Force")={fs};
Physical Surface("Upper")={s1};
Physical Surface("Lower")={s2};
Physical Volume("Column")={V[]};
Physical Line("prova")={prova};






