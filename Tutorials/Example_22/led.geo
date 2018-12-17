L_bar_q = 4; // quantum barrier
Lw = 3;      // QW thickness
Lb = 10;     // barrier thickness


// QW repetition
N = 1;

// well characteristic length
lw = 0.1;
// barrier characteristic length
lb = 0.5;

Point(1) = {-150, 0, 0, 2};
Point(2) = {-L_bar_q, 0, 0, lb};
Point(3) = {0, 0, 0, lw};
Point(4) = {Lw, 0, 0, lw};
x = Lw + L_bar_q;
Point(5) = {x, 0, 0, lb};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Physical Line("buffer") = {1};
Physical Line("well1") = {3};
Physical Line("qbarrier1") = {2,4};

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


x = x + Lb - L_bar_q;
p = newp; Point(p) = {x, 0, 0, 0.5};
l = newl; Line(l) = {p-1, p};
Physical Line("top_barrier") = {l};

// NOTE: last quantum barrier is used as GaN spacer also

x = x + 25;
p = newp; Point(p) = {x, 0, 0, 0.5};
l = newl; Line(l) = {p-1, p};
Physical Line("ebl") = {l};

x = x + 100;
p = newp; Point(p) = {x, 0, 0, 2};
l = newl; Line(l) = {p-1, p};
Physical Line("cap") = {l};

Physical Point("substrate") = {1};
Physical Point("top") = {p};

