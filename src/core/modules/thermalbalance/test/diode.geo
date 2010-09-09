Point(1) = {-5,0,0,0.5};
out[] = Extrude{10,0,0}{ Point{1}; Layers{100};};


Physical Line("Bulk") = {out[1]};
Physical Point("cathode") = {1};
Physical Point("anode") = {out[0]};
