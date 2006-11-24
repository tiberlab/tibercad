#ifndef _WZROTATEDCRYSTAL_H_
#define _WZROTATEDCRYSTAL_H_

#include "RotatedCrystal.h"
class WzRotatedCrystal : public RotatedCrystal
{
 public:

  WzRotatedCrystal();

  WzRotatedCrystal(const double a, const double c);
  
  void set_lat_const(const double a, const double c);


  void calculate_lat_consts();

  void calculate_rot_matrix_miller(std::vector<int> vec_x_mil, std::vector<int> vec_y_mil);

  void set_xyz_mil_direction(std::string dir, int h, int k, int i, int l );


  void read_database (Dummy &db);

 private:

  double a_lat;

  double c_lat;

};


#endif
