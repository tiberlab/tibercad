//Units: nm

//buffer
b = 1e-5;

//Width of InAs layer
wInAs = 3.6;

//Width of GaSb layer
wGaSb = 1.2;

p = newp; Point(p) = {0,0,0,0.1};
Physical Point("left") ={p};
p = newp; Point(p) = {wGaSb,0,0,0.1};
l = newl; Line(l) = {p-1,p};
Physical Line("GaSb") = {l};
p = newp; Point(p) = {wGaSb + wInAs,0,0,0.1};
Physical Point("right") ={p};
l = newl; Line(l) = {p-1,p};
Physical Line("InAs") = {l};

