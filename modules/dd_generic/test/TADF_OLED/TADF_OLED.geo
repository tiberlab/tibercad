

// structrue CYNORA
 
//ITO - aNPD (80 nm) - TCTA (10 nm) - mCBP (10 nm) - [ 80%  mCBP / 20% E001 CT ] (20 nm) - NBPhen (40 nm) - LiF

// anode - HTL - EBL - EM - HBL - ETL - cathode

W_1 = 80;

W_2 = 10;

W_3 = 20;

W_4 = 40;




a1 = 0.1;
a2 =0.25;
a3 = 0.25;

Point(1)=  {0, 0, 0, a3};

anode[] = {1};

t[] = Extrude {W_1, 0, 0} { Point{1}; Layers { W_1/a3 }; Recombine; } ;
HTL[] = t[1]; 

t[] = Extrude {W_2, 0, 0} { Point{news}; Layers { W_2/a2}; Recombine; } ;
EBL[] = t[1]; 

t[] = Extrude {W_3, 0, 0} { Point{news}; Layers { W_3/a2 }; Recombine; } ;
EM[] = t[1]; 

t[] = Extrude {W_4, 0, 0} { Point{news}; Layers { W_4/a3 }; Recombine; } ;
ETL[] = t[1]; 
cathode[] = t[0];


Physical Line("htl") = {HTL[]};
Physical Line("ebl") = {EBL[]};
Physical Line("eml") = {EM[]};
Physical Line("etl") = {ETL[]};
Physical Point("anode") = {anode[]};
Physical Point("cathode") = {cathode[]};




