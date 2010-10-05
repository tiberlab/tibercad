//Length = 1 um

Point(1) = {-1.5,0,0,1};
out[] = Extrude{3,0,0}{ Point{1}; Layers{10};};


Physical Line("Bulk") = {out[1]};
Physical Point("cathode") = {1};
Physical Point("anode") = {out[0]};
