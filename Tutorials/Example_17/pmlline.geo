x = 20;
y = 5;

Point(1) = {0,0,0};

t[] = Extrude {4,0,0} {Point{1}; Layers{4*x}; Recombine; }; 
t[] = Extrude {2,0,0} {Point{2}; Layers{4*x}; Recombine; }; 
t[] = Extrude {2,0,0} {Point{3}; Layers{4*x}; Recombine; }; 
t[] = Extrude {4,0,0} {Point{4}; Layers{4*x}; Recombine; }; 

Physical Line("pml") = {1, 4};
Physical Line("device1") = {2};
Physical Line("device2") = {3};

Physical Point("source") = {2};

Physical Point("left") = {1};
Physical Point("right") = {4};
