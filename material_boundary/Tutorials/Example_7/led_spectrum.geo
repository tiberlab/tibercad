//
//   LED   AlGaN   QW  GaN/AlgAN        AlgaN
//  3 QW  ->   3 separated quantum  regions  






 lc_QW = 0.055 ;    

lc_bigbarr = 0.055;

lc = 10 ; 
lc_b = 1.0  ;


p_algan = 91;  
algan_barr = 12; 
gan_well= 2;  



t1 = 1500; 


t2 = t1 + 25;

start_quantum1 = t2 + (algan_barr/2);  // A

t3 = start_quantum1  + (algan_barr/2) ;



t4 = t3 +gan_well ;

end_quantum1 = t4 + (algan_barr/2);  //  B


t5 = end_quantum1  + (algan_barr/2) ;


t6 = t5 + gan_well;

end_quantum2 = t6 + (algan_barr/2);  // C

t7 = end_quantum2  +  (algan_barr/2)  ;



t8 = t7 +gan_well ;


end_quantum3 = t8 + (algan_barr/2);  //  D

t9 = end_quantum3   +  (algan_barr/2)  ;


t10 = t9 + p_algan;







Point(1) = {0,0,0,lc};


Point(2) = {t1,0,0,lc_bigbarr};


Point(3) = {t2,0,0,lc_QW} ;



Point(4) = {start_quantum1 ,0,0,lc_QW};






Point(5) = {t3,0,0,lc_QW};

Point(6) = { t4 ,0,0,lc_QW};






Point(7) = { end_quantum1 ,0,0,lc_QW};



Point(8) = { t5,0,0,lc_QW};

Point(9) = { t6,0,0,lc_QW};

Point(10) = {end_quantum2 ,0,0,lc_QW};

Point(11) = { t7,0,0,lc_QW};

Point(12) = {t8,0,0,lc_QW};

Point(13) = {end_quantum3,0,0,lc_QW};


Point(14) = {t9,0,0,lc_b};



Point(15) = {t10,0,0,lc_b};


Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};


Line(5) = {5,6};
Line(6) = {6,7};

Line(7) = {7,8};
Line(8) = {8,9};
Line(9) = {9,10};
Line(10) = {10,11};
Line(11) = {11,12};
Line(12) = {12,13};
Line(13) = {13,14};
Line(14) = {14,15};



Physical Point("n_contact") = {1};   // n-contact
Physical Point("p_contact") = {15}; // p-contact


Physical Line("n-algan") = {1};  // n-algan
Physical Line("alinn") = {2}; // alinn
Physical Line("algan_barr1") = {3}; //  Algan barr1


Physical Line("algan_barr_quantum1") = {4,6}; //  Algan  barr  quantum1
Physical Line("GaN_QW1") = {5};   //  GaN QW1

Physical Line("algan_barr_quantum2") = {7,9}; //  Algan  barr quantum2
Physical Line("GaN_QW2") = {8};   //  GaN QW2

Physical Line("algan_barr_quantum3") = {10}; //  Algan  barr quantum3

Physical Line("GaN_QW3") = {11};   //  GaN QW3

Physical Line("i_AlGaN1_quantum") = {12}; //  i-Algan  barr quantum3
Physical Line("i_AlGaN2_noquantum") = {13}; //  i-Algan  barr no quantum

Physical Line("p_AlGaN") = {14};   //  p-AlGaN  




