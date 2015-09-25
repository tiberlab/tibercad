x = 60;

Point(1) = {0,0,0};

t[] = Extrude {300,0,0} {Point{1}; Layers{x}; Recombine; }; 
t[] = Extrude {50,0,0} {Point{2}; Layers{x}; Recombine; }; 
t[] = Extrude {110,0,0} {Point{3}; Layers{x}; Recombine; }; 
t[] = Extrude {35,0,0} {Point{4}; Layers{x}; Recombine; }; 
t[] = Extrude {220,0,0} {Point{5}; Layers{x}; Recombine; }; 
t[] = Extrude {7,0,0} {Point{6}; Layers{x}; Recombine; }; 
t[] = Extrude {200,0,0} {Point{7}; Layers{x}; Recombine; }; 
t[] = Extrude {300,0,0} {Point{8}; Layers{x}; Recombine; }; 

Physical Line("pml") = {1, 8};
Physical Line("air") = {2};
Physical Line("d1") = {3};
Physical Line("d2") = {4};
Physical Line("d3") = {5};
Physical Line("d4") = {6};
Physical Line("d5") = {7};

Physical Point("source") = {2};
