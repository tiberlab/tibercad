cl_channel = 0.6 ; 
cl_ox = 1;
cl_axis = 30; 


l_source = 30; 

l_drain = 30; 

l_intr_s = 20;
l_intr_d = 20;

R_channel = 10; //  NW radius
t_ox = 1.6;


l_channel = 60 ; 

h1= l_channel; 

factor = 1; 


num_lay =   h1/factor;

factor_S = 2; 

factor_D = 2; 

factor_iS =  1;
factor_iD = 1;

num_lay_S = l_source/factor_S; 
num_lay_D = l_source/factor_D; 

num_lay_iS = l_intr_s/factor_iS;
num_lay_iD = l_intr_d/factor_iD;


Point(1) = {0, 0, 0, cl_axis};
Point(2) = {R_channel, 0, 0, cl_channel };
Point(3) = {-R_channel, 0, 0, cl_channel };
Point(4) = {0,R_channel , 0, cl_channel };
Point(5) = {0,-R_channel , 0, cl_channel };



Circle(1) = {2, 1, 4};
Circle(2) = {2, 1, 4};


Circle(3) = {4, 1, 3};
Circle(4) = {3, 1, 5};
Circle(5) = {5, 1, 2};


Point(6) = {R_channel+t_ox, 0, 0, cl_ox };
Point(7) = {-R_channel-t_ox, 0, 0, cl_ox};
Point(8) = {0,R_channel+t_ox , 0, cl_ox };
Point(9) = {0,-R_channel-t_ox , 0, cl_ox };
Circle(6) = {8, 1, 6};
Circle(7) = {8, 1, 7};
Circle(8) = {7, 1, 9};
Circle(9) = {9, 1, 6};


Line(10) = {3, 1};
Line(11) = {1, 2};
Line Loop(12) = {1, 3, 10, 11};
Plane Surface(13) = {12};
Line Loop(14) = {10, 11, -5, -4};
Plane Surface(15) = {14};
Line Loop(16) = {6, -9, -8, -7};
Line Loop(17) = {1, 3, 4, 5};
Plane Surface(18) = {16, 17};





// si  channel a 1

t[] = Extrude   {0.0,0.0,h1/2} { Surface{13}; Layers{ {  20,12 }, { 0.5, 1 } } ; Recombine; } ;

volume = t[1];

Si_volumes[0] = volume;
volume++;


// si  channel a 2


t[] = Extrude   {0.0,0.0,h1/2} { Surface{news - 1}; Layers{ { 12, 20 }, { 0.5, 1 } } ; Recombine; } ;

volume = t[1];

Si_volumes[1] = volume;
volume++;


//drain  intrinsic  region 1 
t[] = Extrude   {0.0,0.0, l_intr_d  } { Surface{news - 1} ; Layers {num_lay_iD}; Recombine; };  

volume = t[1];
drain_intr_volumes[0] = volume;
volume++;

drain_intr_1 = news - 1;



// drain region 1 
t[] = Extrude   {0.0,0.0, l_drain } { Surface{news - 1} ; Layers {num_lay_D}; Recombine; };  
volume = t[1];
drain_volumes[0] = volume;
volume++;

drain_1 = news - 1;



// si  channel b 1

t[] = Extrude   {0.0,0.0,h1/2} { Surface{15}; Layers{ { 20,12 }, { 0.5, 1 }   } ; Recombine; } ;

volume = t[1];
Si_volumes[2] = volume;
volume++;


// si  channel b 2

t[] = Extrude   {0.0,0.0,h1/2} { Surface{news - 1}; Layers{ { 12, 20 }, { 0.5, 1 }   } ; Recombine; } ;

volume = t[1];
Si_volumes[3] = volume;
volume++;



//drain  intrinsic region 2 

t[] = Extrude   {0.0,0.0, l_intr_d  } { Surface{news - 1} ; Layers {num_lay_iD}; Recombine; }; 
 
volume = t[1];
drain_intr_volumes[1] = volume;
volume++;

drain_intr_2 = news - 1;


//drain   region 2 

t[] = Extrude   {0.0,0.0, l_drain } { Surface{news - 1} ; Layers {num_lay_D}; Recombine; };  
volume = t[1];
drain_volumes[1] = volume;
volume++;

drain_2 = news - 1;



//  oxide 1 

t[] = Extrude   {0.0,0.0,h1/2} { Surface{18}; Layers{ { 20,12 }, { 0.5, 1 }   } ; Recombine; } ;

volume = t[1];
Ox_volumes[0] = volume;
volume++;

oxide_surf[0] = t[2]; //
oxide_surf[1] = t[3]; //
oxide_surf[2] = t[4]; //
oxide_surf[3] = t[5]; //



//  oxide 2

t[] = Extrude   {0.0,0.0,h1/2} { Surface{news - 1}; Layers{ { 12,20 }, { 0.5, 1 }   } ; Recombine; } ;

volume = t[1];
Ox_volumes[1] = volume;
volume++;

oxide_surf[4] = t[2]; 
oxide_surf[5] = t[3]; 
oxide_surf[6] = t[4]; 
oxide_surf[7] = t[5]; 



//  source intrinsic region 1

t[] = Extrude   {0.0,0.0,-l_intr_s } { Surface{13}; Layers {num_lay_iS}; Recombine; };  

volume = t[1];
source_intr_volumes[0] = volume;
volume++;

source_intr_1 = news - 1;


//  source  region 1

t[] = Extrude   {0.0,0.0, -l_source  } { Surface{news - 1}; Layers {num_lay_S}; Recombine; };  
volume = t[1];
source_volumes[0] = volume;
volume++;

source_1 = news - 1;


//  source intrinsic region 2

t[] = Extrude   {0.0,0.0,-l_intr_s} { Surface{15}; Layers {num_lay_iS}; Recombine; };  

volume = t[1];
source_intr_volumes[1] = volume;
volume++;

source_intr_2 = news - 1;



//  source region 2

t[] = Extrude   {0.0,0.0,-l_source} { Surface{news - 1}; Layers {num_lay_S}; Recombine; };  
volume = t[1];
source_volumes[1] = volume;
volume++;

source_2 = news - 1;


Physical Volume("channel") = {Si_volumes[], source_intr_volumes[], drain_intr_volumes[] }  ;


Physical Volume("oxide") = Ox_volumes[]  ;


Physical Volume("source_region") = source_volumes[]  ;

Physical Volume("drain_region") = drain_volumes[]  ;



Physical Surface("source") =  {source_1,source_2}; 


Physical Surface("drain") =  {drain_1, drain_2}; 

Physical Surface("gate") = oxide_surf[]; 










