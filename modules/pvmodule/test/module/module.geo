L = 2.5;
H = 0.5;

P1 = 0.02;
P2 = 0.02;
P3 = 0.02;
A = 0.6;

CNT = 0.25;


d = 0.05;

Point(1) = {0, 0, 0, d};
Point(2) = {0, H, 0, d}; 
Point(3) = {L, 0, 0, d};
Point(4) = {L, H, 0, d};
Line(1) = {1, 2};
Line(2) = {3, 4};

x = CNT;
p = newp; Point(p) = {x, 0, 0, d};
p = newp; Point(p) = {x, H, 0, d};
l = newl; Line(l) = {1, p-1};
l = newl; Line(l) = {p-1, p};
l = newl; Line(l) = {p, 2};
ll = news; Line Loop(ll) = {l-2, l-1, l, -1};
Plane Surface(ll) = {ll};
p3_reg[] = {ll};

x = CNT + P2;
p = newp; Point(p) = {x, 0, 0, d};
p = newp; Point(p) = {x, H, 0, d};
l = newl; Line(l) = {p-3, p-1};
l = newl; Line(l) = {p-1, p};
l = newl; Line(l) = {p, p-2};
ll = news; Line Loop(ll) = {l-2, l-1, l, 5-l};
Plane Surface(ll) = {ll};
p2_reg[] = {ll};

x = x + P1;
p = newp; Point(p) = {x, 0, 0, d};
p = newp; Point(p) = {x, H, 0, d};
l = newl; Line(l) = {p-3, p-1};
l = newl; Line(l) = {p-1, p};
l = newl; Line(l) = {p, p-2};
ll = news; Line Loop(ll) = {l-2, l-1, l, 5-l};
Plane Surface(ll) = {ll};
p1_reg[] = {ll};

x = x + A;
p = newp; Point(p) = {x, 0, 0, d};
p = newp; Point(p) = {x, H, 0, d};
l = newl; Line(l) = {p-3, p-1};
l = newl; Line(l) = {p-1, p};
l = newl; Line(l) = {p, p-2};
ll = news; Line Loop(ll) = {l-2, l-1, l, 5-l};
Plane Surface(ll) = {ll};
act_reg[] = {ll};

// second stripe

x = x + P3;
p = newp; Point(p) = {x, 0, 0, d};
p = newp; Point(p) = {x, H, 0, d};
l = newl; Line(l) = {p-3, p-1};
l = newl; Line(l) = {p-1, p};
l = newl; Line(l) = {p, p-2};
ll = news; Line Loop(ll) = {l-2, l-1, l, 5-l};
Plane Surface(ll) = {ll};
p3_reg[] = {p3_reg[], ll};

x = x + P2;
p = newp; Point(p) = {x, 0, 0, d};
p = newp; Point(p) = {x, H, 0, d};
l = newl; Line(l) = {p-3, p-1};
l = newl; Line(l) = {p-1, p};
l = newl; Line(l) = {p, p-2};
ll = news; Line Loop(ll) = {l-2, l-1, l, 5-l};
Plane Surface(ll) = {ll};
p2_reg[] = {p2_reg[], ll};

x = x + P1;
p = newp; Point(p) = {x, 0, 0, d};
p = newp; Point(p) = {x, H, 0, d};
l = newl; Line(l) = {p-3, p-1};
l = newl; Line(l) = {p-1, p};
l = newl; Line(l) = {p, p-2};
ll = news; Line Loop(ll) = {l-2, l-1, l, 5-l};
Plane Surface(ll) = {ll};
p1_reg[] = {p1_reg[], ll};

x = x + A;
p = newp; Point(p) = {x, 0, 0, d};
p = newp; Point(p) = {x, H, 0, d};
l = newl; Line(l) = {p-3, p-1};
l = newl; Line(l) = {p-1, p};
l = newl; Line(l) = {p, p-2};
ll = news; Line Loop(ll) = {l-2, l-1, l, 5-l};
Plane Surface(ll) = {ll};
act_reg[] = {act_reg[], ll};

// second stripe

x = x + P3;
p = newp; Point(p) = {x, 0, 0, d};
p = newp; Point(p) = {x, H, 0, d};
l = newl; Line(l) = {p-3, p-1};
l = newl; Line(l) = {p-1, p};
l = newl; Line(l) = {p, p-2};
ll = news; Line Loop(ll) = {l-2, l-1, l, 5-l};
Plane Surface(ll) = {ll};
p3_reg[] = {p3_reg[], ll};

x = x + P2;
p = newp; Point(p) = {x, 0, 0, d};
p = newp; Point(p) = {x, H, 0, d};
l = newl; Line(l) = {p-3, p-1};
l = newl; Line(l) = {p-1, p};
l = newl; Line(l) = {p, p-2};
ll = news; Line Loop(ll) = {l-2, l-1, l, 5-l};
Plane Surface(ll) = {ll};
p2_reg[] = {p2_reg[], ll};

x = x + P1;
p = newp; Point(p) = {x, 0, 0, d};
p = newp; Point(p) = {x, H, 0, d};
l = newl; Line(l) = {p-3, p-1};
l = newl; Line(l) = {p-1, p};
l = newl; Line(l) = {p, p-2};
ll = news; Line Loop(ll) = {l-2, l-1, l, 5-l};
Plane Surface(ll) = {ll};
p1_reg[] = {p1_reg[], ll};

x = x + A;
p = newp; Point(p) = {x, 0, 0, d};
p = newp; Point(p) = {x, H, 0, d};
l = newl; Line(l) = {p-3, p-1};
l = newl; Line(l) = {p-1, p};
l = newl; Line(l) = {p, p-2};
ll = news; Line Loop(ll) = {l-2, l-1, l, 5-l};
Plane Surface(ll) = {ll};
act_reg[] = {act_reg[], ll};


l = newl; Line(l) = {p-1, 3};
l = newl; Line(l) = {4, p};
ll = news; Line Loop(ll) = {l-1, 2, l, 4-l};
Plane Surface(ll) = {ll};
p3_reg[] = {p3_reg[], ll};



Physical Surface("active") = {act_reg[]};
Physical Surface("P1") = {p1_reg[]};
Physical Surface("P2") = {p2_reg[]};
Physical Surface("P3") = {p3_reg[]};
Physical Line("anode") = {1};
Physical Line("cathode") = {2};

//N = 11;
//Transfinite Curve {1} = N Using Progression 1;
//Transfinite Curve {2} = N Using Progression 1;
//Transfinite Curve {3} = N Using Progression 1;
//Transfinite Curve {4} = N Using Progression 1;
//Transfinite Surface {1};
//Recombine Surface {1};
