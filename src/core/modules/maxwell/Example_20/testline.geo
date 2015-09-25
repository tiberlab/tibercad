x = 60;
l = 3.14159265 * 299792458;

Point(1) = {0,0,0};

t[] = Extrude {2*l,0,0} {Point{1}; Layers{x}; Recombine; }; 
t[] = Extrude {l,0,0} {Point{2}; Layers{x}; Recombine; }; 
t[] = Extrude {2*l,0,0} {Point{3}; Layers{x}; Recombine; }; 

Physical Line("device") = {2};

Physical Line("pml") = {1, 3};

Physical Point("start") = {1};
Physical Point("end") = {4};
