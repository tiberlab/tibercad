cGaN = 0.5185;
cAlN = 0.4982;

NGaN = 2;
NAlN = 1;

LGaN = cGaN*NGaN + 1e-6;
LAlN = cGaN*NAlN + 1e-6;

N = 8;
x = 0;

aln[] = {};
gan[] = {};

p = newp; Point(p) = {0,0,0,0.1};

For i In {0:2*N}
 x = x + LGaN;
 p = newp; Point(p) = {x, 0, 0, 0.1};
 x = x + LAlN;
 p = newp; Point(p) = {x, 0, 0, 0.1};

 l = newl; Line(l) = {p-2, p-1};
 l = newl; Line(l) = {p-1, p};

 If (i == N)
  Physical Line("qGaN") = {l-1};
  Physical Line("qAlN") = {l};
 Else
  gan[] = {gan[], l-1};
  aln[] = {aln[], l};
 EndIf

EndFor

Physical Line("GaN") = {gan[]};
Physical Line("AlN") = {aln[]};
Physical Point("left") ={1};
Physical Point("right") ={p};

