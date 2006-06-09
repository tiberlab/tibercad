//width = 20.0;
width = 40.0;
well = 3.0;
length = 100 + well;
lc = width / 4.0;

steps_doped = 1;
steps_well = 50;

xmax = width / 2;
xmin = -xmax;
ymax = length / 2;
ymin = -ymax;


Point(1) = {xmin, ymin, 0, lc};
Point(2) = {xmax, ymin, 0, lc};
s = news; Line(s) = {1,2};

N = 200;

l_n = 0.5 - 0.5 * well / length;
vol = 1;
For i In {0:N-2}
  v_n_gan[i] = vol;
  height[i] = l_n * Log10(i+2) / Log10(N);
  steps[i] = steps_doped;
  layers[i] = vol;
  vol = vol + 1;
EndFor


v_well = vol;
height[vol-1] = height[vol-2] + well / length;
steps[vol-1] = steps_well;
layers[vol-1] = vol;
vol = vol + 1;

N0 = vol - 1;
For i In {0:N-2}
  v_p_gan[i] = vol;
  height[i+N0] = 1 - l_n * Log10(N-i-1) / Log10(N);
  steps[i+N0] = steps_doped;
  layers[i+N0] = vol;
  vol = vol + 1;
EndFor


Extrude Line {s, {0.0, length, 0.0}} {Recombine; Layers{steps[], layers[], height[]};};

Physical Surface(1) = {v_n_gan[]};
Physical Surface(2) = {v_well};
Physical Surface(3) = {v_p_gan[]};
Physical Line(1) = {2};
Physical Line(2) = {1};
