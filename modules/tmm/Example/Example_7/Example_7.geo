L_Air = 1000;

Point(1) = {0,0,0,1};
Point(2) = {L_Air,0,0,1};

Line(1) = {1,2};

Physical Line("Air") = {1};










