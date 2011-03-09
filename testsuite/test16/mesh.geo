SideX = 10;
SideZ = 10;
Dx = 53/3;
Dz = 53/3;
Nx = 1;
Nz = 2;
H = 10;

W = Nx  * Dx;
L = Nz * Dz;

clb = 3;
clp =1;

Point(0) = {0, 0, 0, clb};
Point(1) = {W, 0, 0, clb};
Point(2) = {0, H, 0, clb};
Point(3) = {W, H, 0, clb};
Point(4) = {0, 0, L, clb};
Point(5) = {W, 0, L, clb};
Point(6) = {0, H, L, clb};
Point(7) = {W, H, L, clb};
Line(1) = {6, 4};
Line(2) = {4, 5};
Line(3) = {5, 7};
Line(4) = {7, 6};
Line(5) = {6, 2};
Line(6) = {2, 0};   
Line(7) = {0, 1};
Line(8) = {1, 3};
Line(9) = {3, 2};
Line(10) = {1, 5};
Line(11) = {7, 3};
Line(12) = {0, 4};
Line Loop(14) = {4, 1, 2, 3};
Plane Surface(14) = {14};
Line Loop(16) = {9, 6, 7, 8};
Plane Surface(16) = {16};
Line Loop(18) = {11, -8, 10, 3};
Plane Surface(18) = {18};
Line Loop(20) = {12, -1, 5, 6};
Plane Surface(20) = {20};

//Initial Line Loop
UpLL[0] = newl; Line Loop(UpLL[0]) = {-5,-4,11,9};
LoLL[0] = newl; Line Loop(LoLL[0]) = {12,2,-10,-7};
Nup = 0;
Nlo = 0;
//Lateral surfaces
S[0] = 14;
S[1] = 16;
S[2] = 18;
S[3] = 20;
Ns = 3;
//----------------
Np = 0;



For tx In {1:Nx}

PosX = Dx/2 + (tx-1)*Dx;

 For tz In {1:Nz}

    PosZ = Dz/2 + (tz-1)*Dz;


//--------Make a square pore--------

p0=newp; Point(p0) = {PosX-SideX/2, 0,PosZ-SideZ/2, clp};
p1=newp; Point(p1) = {PosX+SideX/2, 0,PosZ-SideZ/2, clp};
p2=newp; Point(p2) = {PosX+SideX/2, 0,PosZ+SideZ/2, clp};
p3=newp; Point(p3) = {PosX-SideX/2, 0,PosZ+SideZ/2, clp};
p4=newp; Point(p4) = {PosX-SideX/2, H,PosZ-SideZ/2, clp};
p5=newp; Point(p5) = {PosX+SideX/2, H,PosZ-SideZ/2, clp};
p6=newp; Point(p6) = {PosX+SideX/2, H,PosZ+SideZ/2, clp};
p7=newp; Point(p7) = {PosX-SideX/2, H,PosZ+SideZ/2, clp};
l0 = newl; Line(l0) = {p0,p1};
l1 = newl; Line(l1) = {p1,p2};
l2 = newl; Line(l2) = {p2,p3};
l3 = newl; Line(l3) = {p3,p0};
l4 = newl; Line(l4) = {p4,p5};
l5 = newl; Line(l5) = {p5,p6};
l6 = newl; Line(l6) = {p6,p7};
l7 = newl; Line(l7) = {p7,p4};
l8 = newl; Line(l8) = {p0,p4};
l9 = newl; Line(l9) = {p5,p1};
l10 = newl; Line(l10) = {p3,p7};
l11 = newl; Line(l11) = {p6,p2};

//Upper Line Loop
Nup = Nup + 1;
UpLL[Nup]=newl; Line Loop(UpLL[Nup]) = {l4, l5, l6, l7};
//CubeS[0]=news; Plane Surface(CubeS[0]) = {UpLL[Nup]};

//Lower Line Loop
Nlo = Nlo + 1;
LoLL[Nlo]=newl; Line Loop(LoLL[Nlo]) = {l0, l1, l2, l3};
//CubeS[1]=news; Plane Surface(CubeS[1]) = {LoLL[Nlo]};

//Lateral pore surfaces
ll=newl; Line Loop(ll) = {l2, l10, -l6, l11};
Ns = Ns + 1;
S[Ns] = news;
CubeS[2] = S[Ns];
Plane Surface(S[Ns]) = {ll};

//Lateral pore surfaces
ll=newl; Line Loop(ll) = {l0, -l9, -l4, -l8};
Ns = Ns + 1;
S[Ns] = news;
CubeS[3] = S[Ns];
Plane Surface(S[Ns]) = {ll};

//Lateral pore surfaces
ll=newl; Line Loop(ll) = {l1, -l11, -l5, l9};
Ns = Ns + 1;
S[Ns] = news;
CubeS[4] = S[Ns];
Plane Surface(S[Ns]) = {ll};

//Lateral pore surfaces
ll=newl; Line Loop(ll) = {-l3, l10,l7,-l8};
Ns = Ns + 1;
S[Ns] = news;
CubeS[5] = S[Ns];
Plane Surface(S[Ns]) = {ll};

//Volume of the pore
//SV = news; Surface Loop(SV) = {CubeS[]};
//[Np] = newv; Volume(V[Np]) = {SV};
Np = Np + 1;
//--------------------------------------

EndFor

EndFor



//Create upper surface
Ns = Ns + 1;
S[Ns] = news;
Plane Surface(S[Ns]) = {UpLL[]};

//Create lower surface
Ns = Ns + 1;
S[Ns] = news;
Plane Surface(S[Ns]) = {LoLL[]};

//Create the Bulk Volume
SL = news; Surface Loop(SL) = {S[]};
VB = newv; Volume(VB) = {SL};

//-----
Physical Volume("Bulk") = {VB};
//Physical Volume("Pore") = {V[]};
Physical Surface("Hot") = {S[0]};
Physical Surface("Cold") = {S[1]};
Physical Surface("Lower") = {S[Ns]};





