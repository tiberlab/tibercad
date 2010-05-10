//Units: A

well_lenght = 50.0;

quantum_barrier_lenght = 200.0;

device_lenght = 500.0;

lc_QW = 0.2;
lc_barrier = 0.4;
lc_device = 1.0;


Point(1) = {0,0,0,lc_device};
Point(2) = {device_lenght,0, 0, lc_device};

Point(3) = {device_lenght + quantum_barrier_lenght, 0, 0, lc_barrier};
Point(4) = {device_lenght + quantum_barrier_lenght + well_lenght, 0, 0, lc_QW};

Point(5) = {device_lenght + 2*quantum_barrier_lenght + well_lenght,0,0,lc_barrier};
Point(6) = {2*device_lenght + 2*quantum_barrier_lenght + well_lenght,0,0,lc_device};


Physical Point("cathode") = {1};
Physical Point("anode") = {6};

Line(1) = {3,4};
Physical Line("QW") = {1};

Line(2) = {1,2};
Physical Line("device_left") = {2};

Line(3) = {2,3};
Physical Line("barrier_left") = {3};

Line(4) = {4,5};
Physical Line("barrier_right") = {4};

Line(5) = {5,6};
Physical Line("device_right") = {5};
