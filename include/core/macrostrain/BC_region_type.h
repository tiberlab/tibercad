class BC_region_type
{
 public:
  Real coord_min[3];
  Real coord_max[3];
  int  material; 
};



class stress_region_type
{
 public:
  Real coord_min[3];
  Real coord_max[3];
  double  stress; 
};
