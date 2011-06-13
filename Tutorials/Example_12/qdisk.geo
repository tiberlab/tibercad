s = 0.31736; // hexagon side length

a = 63.5 * s;  // 7

d = a / 9.2; // 1

d1 = 1;
hw = 3;
hb = 4; // 5
h1 = 80; // 100
h2 = 80; // 100
h3=50;
N=4;

Function add_GaN_side_surfaces

   interface_plane = news - 1;
   For i1 In {0:5}
       GaN_side_surface_number++;
       GaN_side_surface[GaN_side_surface_number] = interface_plane - 1 - i1*4;
   EndFor

Return


Function add_AlGaN_side_surfaces

   interface_plane = news - 1;
   For i2 In {0:5}
       AlGaN_side_surface_number++;
       AlGaN_side_surface[AlGaN_side_surface_number] = interface_plane - 1 - i2*4;
   EndFor

Return


Point(1) = {Sqrt(3)/2*a, a/2, 0, d};
Point(2) = {0,a,0,d};
Point(3) = {-Sqrt(3)/2*a,a/2,0,d};
Point(4) = {-Sqrt(3)/2*a,-a/2,0,d};
Point(5) = {0,-a,0,d};
Point(6) = {Sqrt(3)/2*a, -a/2, 0, d};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};
Line(6) = {6,1};

Line Loop(7) = {1,2,3,4,5,6};
Plane Surface(1) = {7};


N = 30;
For i In {0:N-2}
  heightb[i] = 0.5 * (1 + Log10(i+2) / Log10(N));
  stepsb[i] = 1;
  heightt[i] = 0.5 * (1 - Log10(N-i-1) / Log10(N));
  stepst[i] = 1;
EndFor
heightt[N-2] = heightt[N-2] - 0.05;



//-----------------------------
// GaN
GaN_side_surface_number =  0;
AlGaN_side_surface_number =  0;
volume = 1;
//------------n_part------------------------------------------------------------------






t[] = Extrude   {0.0,0.0,h1} { Surface{1}; Layers { { 10, stepsb[] }, { 0.55, heightb[] } }; Recombine; };  //AlGaN 
volume = t[1];
AlGaN_volumes[0] = volume;
Call add_AlGaN_side_surfaces;
volume++;






t[] = Extrude {0.0,0.0,hb} {Surface{news - 1}; Layers {12}; Recombine;} ;//AlGaN 
volume = t[1];
AlGaN_volumes[1] = volume;
Call add_AlGaN_side_surfaces;
volume++;

//----i_part--------------------------------------------------------------------------



t[] = Extrude {0.0,0.0,hw} {Surface{news - 1};Layers{10}; Recombine;};//GaN
volume = t[1];
GaN_volumes[0] =  volume;
Call add_GaN_side_surfaces;
volume++;


//----p_part------------------------------------------------------------------------



t[] = Extrude {0.0,0.0,hb} {Surface{news - 1}; Layers {12}; Recombine;} ;//AlGaN 
volume = t[1];

AlGaN_volumes[2] = volume;
Call add_AlGaN_side_surfaces;
volume++;



t[] = Extrude {0.0,0.0,h1} {Surface{news - 1}; Layers { { stepst[], 10 }, { heightt[], 1 }  }; Recombine;} ;//AlGaN 
volume = t[1];
AlGaN_volumes[3] = volume;
Call add_AlGaN_side_surfaces;
volume++;





//-----------end_of_buffer-----------------------------------------------------------

//----------------------------------------------------------------------------------------
//Physical volumes------------------------------------------------------------------------

Physical Volume("well") =GaN_volumes[0] ;
Physical Volume("n-side") =AlGaN_volumes[0];
Physical Volume("barrier1") =AlGaN_volumes[1];
Physical Volume("barrier2") =AlGaN_volumes[2];
Physical Volume("p-side") =AlGaN_volumes[3];








//------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
//Physical surfaces------------------------------------------------------------------------

Physical Surface("anode") =  {news - 1}; //p_contact
Physical Surface("cathode") = {1}; //n_contact
Physical Surface("GaN_surface") =  GaN_side_surface[];    //GaN side surfaces
Physical Surface("AlGaN_surface") =  AlGaN_side_surface[];  //AlGaN side surfaces




