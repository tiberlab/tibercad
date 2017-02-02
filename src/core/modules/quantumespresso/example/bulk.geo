x = 0.0;
c = x * 0.5703 + (1 - x) * 0.5185;

N = 2;
L = c*(N-0.1);
Point(1) = {0,0,0,0.01};
Point(2) = {L, 0, 0, 0.01};

Line(1) = {1, 2};

Physical Line("bulk") ={1};
Physical Point("left") ={1};
Physical Point("right") ={2};

