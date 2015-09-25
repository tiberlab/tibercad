x = 60;
l = 3.14159265 * 299792458;

Point(1) = {0,0,0};

t[] = Extrude {l,0,0} {Point{1}; Layers{x}; Recombine; }; 

Physical Line("air") = {1};

Physical Point("start") = {1};
Physical Point("end") = {2};
