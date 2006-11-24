#ifndef _ZBROTATEDCRYSTAL_H_
#define _ZBROTATEDCRYSTAL_H_

#include "RotatedCrystal.h"

class ZbRotatedCrystal : public RotatedCrystal
{

 public:
  ZbRotatedCrystal();

  ZbRotatedCrystal(const double a);

  void set_lat_const(const double a);

  virtual void calculate_lat_consts();

  void calculate_rot_matrix_miller(std::vector<int> vec_x_mil, std::vector<int> vec_y_mil);

  void set_xyz_mil_direction(std::string dir, int h, int k, int l);

  virtual void read_database (Dummy &db);

 

 private:

  //! Lattice constant
  double a_lat; 


};

#endif
