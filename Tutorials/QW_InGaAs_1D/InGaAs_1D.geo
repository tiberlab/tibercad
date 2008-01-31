//
//  QW  GaAs/InGaAs/GaAs 



lc_QW = 0.20;  



lc = 1.5 ;   



end_gan_buffer = 500  ; 
first_barrier_thick = 100 ;  

end_first_barrier = end_gan_buffer + first_barrier_thick  ;

start_quantum_region = end_first_barrier -5 ;  //

well_width = 8.0; 
  //4.00  ;


end_QW = end_first_barrier + well_width;

end_quantum_region = end_QW  +  5 ;  


top_barrier_thick = 10;  
end = end_QW + top_barrier_thick ;

Point(1) = {0,0,0,lc};

Point(2) = {end_gan_buffer,0,0,lc};

Point(3) = {start_quantum_region,0,0,lc_QW} ;

Point(4) = {end_first_barrier,0,0,lc_QW};
Point(5) = {end_QW   ,0,0,lc_QW};

Point(6) = { end_quantum_region  ,0,0,lc_QW};

Point(7) = {end ,0,0,lc};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};
Line(6) = {6,7};


Physical Point(1) = {1};
Physical Point(2) = {7};

Physical Line(1) = {1};
Physical Line(2) = {2};
Physical Line(3) = {3};
Physical Line(4) = {4};
Physical Line(5) = {5};
Physical Line(6) = {6};


