SetFactory("OpenCASCADE");


W = 25;
Hb = 25;
Hd = 2;
Hc = 0.7;
Lg = 5;
Ld = 5;
d = 1;
dc = 0.2;


Rectangle(1) = {-W/2, -Hb, 0, W, Hb, 0};
Rectangle(2) = {-W/2, 0, 0, W, Hc, 0};
Rectangle(3) = {-Lg/2, Hc, 0, Lg, Hd, 0};

p = newp; Point(p) = {-W/2 + Ld, Hc, 0, dc};
p = newp; Point(p) = {W/2 - Ld, Hc, 0, dc};


BooleanFragments{ Surface{2}; Delete; }{ Point{p, p-1}; Delete; }
BooleanFragments{ Surface{1}; Delete; }{ Surface{2, 3}; Delete; }

l() = Unique(Abs(Boundary{ Surface{2}; }));
//Printf("ids: ", l());
Physical Line("source") = {l(3)};
Physical Line("drain") = {l(4)};

p() = Unique(Abs(Boundary{ Line{l(3), l(4), l(5), l(6), l(7)}; }));
Characteristic Length {p()} = dc;
p() = Unique(Abs(Boundary{ Line{l(0)}; }));
Characteristic Length {p()} = d;

l() = Unique(Abs(Boundary{ Surface{3}; }));
Physical Line("gate") = {l(3)};


Physical Surface("substrate") = {1};
Physical Surface("channel") = {2};
Physical Surface("dielectric") = {3};
Physical Line("back") = {1};
