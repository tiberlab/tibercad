dgan=7.0;
dqw=3.0;
dn=100;
dp=25;

Point(1) = {0, 0, 0, 0.5};
Point(2) = {dn, 0, 0, 0.2};
Point(3) = {dn+0*dqw+1*dgan, 0, 0, 0.05};
Point(4) = {dn+1*dqw+1*dgan, 0, 0, 0.05};
Point(5) = {dn+1*dqw+2*dgan, 0, 0, 0.05};
Point(6) = {dn+2*dqw+2*dgan, 0, 0, 0.05};

Point(7) = {dn+2*dqw+2.5*dgan, 0, 0, 0.05};
Point(8) = {dn+2*dqw+3*dgan, 0, 0, 0.05};
Point(9) = {dn+3*dqw+3*dgan, 0, 0, 0.05};
Point(10) = {dn+3*dqw+3.5*dgan, 0, 0, 0.05};

Point(11) = {dn+3*dqw+4*dgan, 0, 0, 0.05};
Point(12) = {dn+4*dqw+4*dgan, 0, 0, 0.05};
Point(13) = {dn+4*dqw+5*dgan, 0, 0, 0.05};
Point(14) = {dn+5*dqw+5*dgan, 0, 0, 0.05};
Point(15) = {dn+5*dqw+6*dgan, 0, 0, 0.05};
Point(16) = {dn+5*dqw+6*dgan+dp, 0, 0, 0.5};

Line (1) = {1, 2};
Line (2) = {2, 3};
Line (3) = {3, 4};
Line (4) = {4, 5};
Line (5) = {5, 6};
Line (6) = {6, 7};
Line (7) = {7, 8};
Line (8) = {8, 9};
Line (9) = {9, 10};
Line (10) = {10, 11};
Line (11) = {11, 12};
Line (12) = {12, 13};
Line (13) = {13, 14};
Line (14) = {14, 15};
Line (15) = {15, 16};


Physical Line  ("nGaN") = {1};
Physical Line  ("well") = {3,5,11,13};
Physical Line  ("well_q") = {8};
Physical Line  ("barrier") = {2,4,6,10,12,14};
Physical Line  ("barrier_q") = {7,9};
Physical Line  ("pGaN") = {15};

Physical Point ("substrate") = {1};
Physical Point ("top") = {16};
