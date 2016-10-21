L = 2;
H = 1;
d = 0.03;
dc = 50;

Point(1) = {0, 0, 0, d};
Point(2) = {L, 0, 0, d};

Line(1) = {1, 2};

t[] = Extrude {0, H, 0} {Line{1}; Layers{dc};};

Physical Line('anode') = {2};
Physical Line('cathode') = {1};
Physical Surface('bulk') = {5};
Physical Line('no_contact') = {3};

