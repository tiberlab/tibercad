scale = 1000;

 la=0.02 * scale;
 dim = 1 * scale;
Point(1) = {0,0,0,la};
Point(2) = {dim/2,0,0,la};
Point(3) = {dim,0,0,la};

Line(1) = {1,2};
Line(2) = {2,3};
 

  Physical Line(1) = {1};
  Physical Line(2) = {2};
  Physical Point(1) = {1};
  Physical Point(2) = {3};
  Physical Point(3) = {2};
    
    
  

