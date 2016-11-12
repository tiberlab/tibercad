L_bar_q = 4; // quantum barrier
Lw = 3;      // QW thickness
Lb = 20;     // barrier thickness


// number of QWs
N = 1;

// well characteristic length
lw = 0.1;
// barrier characteristic length
lb = 0.2;

Point(1) = {-150, 0, 0, 2};
Point(2) = {-2*L_bar_q, 0, 0, lb};
Point(3) = {-L_bar_q, 0, 0, lb};
Point(4) = {0, 0, 0, lw};
Point(5) = {Lw, 0, 0, lw};
x = Lw + L_bar_q;
Point(6) = {x, 0, 0, lb};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};
Physical Line("buffer") = {1};
Physical Line("buffer_q") = {2};
Physical Line("well1") = {4};
Physical Line("qbarrierl") = {3};
Physical Line("qbarriert") = {5};

For i In {2:N}
  x = x + Lb - 2*L_bar_q;
  p = newp; Point(p) = {x, 0, 0, lb};
  l = newl; Line(l) = {p-1, p};
  Physical Line(Sprintf("barrier%g", i-1)) = {l};
  x = x + L_bar_q;
  p = newp; Point(p) = {x, 0, 0, lw};
  l = newl; Line(l) = {p-1, p};
  x = x + Lw;
  p = newp; Point(p) = {x, 0, 0, lw};
  l = newl; Line(l) = {p-1, p};
  Physical Line(Sprintf("well%g", i)) = {l};
  x = x + L_bar_q;
  p = newp; Point(p) = {x, 0, 0, lb};
  l = newl; Line(l) = {p-1, p};
  Physical Line(Sprintf("qbarrier%g", i)) = {l-2,l};
EndFor



// NOTE: last quantum barrier is used as GaN spacer also

x = x + 30;
p = newp; Point(p) = {x, 0, 0, 0.5};
l = newl; Line(l) = {p-1, p};
Physical Line("ebl") = {l};

x = x + 50;
p = newp; Point(p) = {x, 0, 0, 2};
l = newl; Line(l) = {p-1, p};

Physical Line("cap") = {l};

x = x + 50;
p = newp; Point(p) = {x, 0, 0, 2};
l = newl; Line(l) = {p-1, p};

Physical Line("p+") = {l};

Physical Point("substrate") = {1};
Physical Point("top") = {p};


