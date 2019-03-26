l = 0.6;
d = 0.1;


Point(1) = {0, 0, 0, d};
Point(2) = {l, 0, 0, d};


Line(1) = {1,2};


Physical Line("bulk") = {1};

