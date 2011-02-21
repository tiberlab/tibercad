scale = 30;
wb=5*scale;
hb=2*scale;
wt=wb/2;
ht=1*scale;
ct=1*scale;
air_x=5*scale;
air_y=7*scale;
sub_x = -2 * scale;



ms  =2;

c_cont1 = 0.001 * scale * ms;
c_cont2 = 0.02 * scale * ms;
c_jun = 0.03 * scale * ms;
c_sub = 0.2 * scale * ms;
c_sub2 = 0.05 * scale * ms ;
c_sub_long = 2 *scale * ms ;
c_air = 1 * scale * ms;


p1=newp;Point(p1) = {0, 0, 0, c_sub};
p2=newp;Point(p2) = {wb, 0, 0, c_sub};
p3=newp;Point(p3) = {wb, hb, 0, c_cont2};
p4=newp;Point(p4) = {wb-ct, hb, 0, c_cont2};
p5=newp;Point(p5) = {wt, hb, 0, c_jun};
p6=newp;Point(p6) = {0, hb, 0, c_jun};
p7=newp;Point(p7) = {wt, hb+ht, 0, c_sub2};
p8=newp;Point(p8) = {ct, hb+ht, 0, c_cont1};
p9=newp;Point(p9) = {0, hb+ht, 0, c_cont1};
p10=newp;Point(p10) = {wb+air_x,0, 0, c_air};
p11=newp;Point(p11) = {0, hb+ht+air_y, 0, c_air};
p12=newp;Point(p12) = {wb+air_x,hb+ht+air_y,0,c_air};
p13=newp;Point(p13) = {0, hb+ht*4/5, 0, c_sub2};
p14=newp;Point(p14) = {wt,hb+ht*4/5,0,c_sub2};





l1=newl;Line(l1) = {p1,p2};
l2=newl;Line(l2) = {p2,p3};
l3=newl;Line(l3) = {p3,p4};
l4=newl;Line(l4) = {p4,p5};
l5=newl;Line(l5) = {p5,p6};
l6=newl;Line(l6) = {p6,p1};
l7=newl;Line(l7) = {p5,p14};
l8=newl;Line(l8) = {p14,p7};
l9=newl;Line(l9) = {p7,p8};
l10=newl;Line(l10) = {p8,p9};
l11=newl;Line(l11) = {p9,p13};
l12=newl;Line(l12) = {p13,p6};
l13=newl;Line(l13) = {p13,p14};
l14=newl;Line(l14) = {p9,p11};
l15=newl;Line(l15) = {p11,p12};
l16=newl;Line(l16) = {p12,p10};
l17=newl;Line(l17) = {p10,p2};



l23=newl;Line Loop(l23) = {l1,l2,l3,l4,l5,l6};
s1=news;Plane Surface(s1) = {l23};

l24=newl;Line Loop(l24) = {l5,-l12,l13,-l7};
s2=news;Plane Surface(s2) = {l24};

l25=newl;Line Loop(l25) = {l13,l8,l9,l10,l11};
s3=news;Plane Surface(s3) = {l25};

l26=newl;Line Loop(l26) = {l14,l15,l16,l17,l2,l3,l4,l7,l8,l9,l10};
s4=news;Plane Surface(s4) = {l26};


Physical Surface("n_side") = {s1}; // nside
Physical Surface("p_side") = {s2,s3}; // pside
Physical Surface("air") = {s4}; // thermal



lp1=newl;Physical Line("anode") = {l10}; // anode
lp2=newl;Physical Line("cathode") = {l3}; // cathode
Physical Line("reservoir") = {l15,l16,l17}; // thermal
Physical Line("dissipator") = {l1}; // thermal
