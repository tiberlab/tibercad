x = 60;

Point(1) = {0,0,0};

t[] = Extrude {300,0,0} {Point{1}; Layers{x}; Recombine; }; 
t[] = Extrude {50,0,0} {Point{2}; Layers{x}; Recombine; }; 
t[] = Extrude {100,0,0} {Point{3}; Layers{x}; Recombine; }; 
t[] = Extrude {150,0,0} {Point{4}; Layers{x}; Recombine; }; 
t[] = Extrude {200,0,0} {Point{5}; Layers{x}; Recombine; }; 
t[] = Extrude {280,0,0} {Point{6}; Layers{x}; Recombine; }; 
t[] = Extrude {5,0,0} {Point{7}; Layers{x}; Recombine; }; 
t[] = Extrude {200,0,0} {Point{8}; Layers{x}; Recombine; }; 

Physical Line("pml") = {1};
Physical Line("air") = {2};
Physical Line("sio2") = {3};
Physical Line("ito") = {4};
Physical Line("pedot") = {5};
Physical Line("blend") = {6};
Physical Line("Ca") = {7};
Physical Line("Al") = {8};

Physical Point("source") = {2};
