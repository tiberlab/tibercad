


start_qw = -10;  
end_qw = 10; 


Point(1) = {-100, 0, 0, 4};
Point(2) = {-20, 0, 0, 1.0};
Point(3) = {start_qw, 0, 0, 0.1};
Point(4) = {end_qw, 0, 0, 0.1};
Point(5) = {20, 0, 0, 1.0};
Point(6) = {100, 0, 0, 4};



Line (1) = {1, 2};
Line (2) = {2, 3};
Line (3) = {3, 4};
Line (4) = {4, 5};
Line (5) = {5, 6};


Physical Line ("well") = {3};
Physical Line ("barrier1_a") = {1};
Physical Line ("barrier2_a") = {5};
Physical Line ("barrier1_q") = {2};
Physical Line ("barrier2_q") = {4};



Physical Point ("substrate") = {1};
Physical Point ("top") = {6};
