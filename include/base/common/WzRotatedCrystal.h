#ifndef _WZROTATEDCRYSTAL_H_
#define _WZROTATEDCRYSTAL_H_

#include "RotatedCrystal.h"
class WzRotatedCrystal : public RotatedCrystal
{
 public:
  //!Constructor
  /*!
    sets growth directions
    x [ 1, 0, -1, 0]
    y [-1, 2, -1, 0]
    z [ 0, 0,  0, 1]
   */
  WzRotatedCrystal(const ModelOptions& options);

  //WzRotatedCrystal(const double a, const double c);
  
  void set_lat_const(const double a, const double c);


  void calculate_lat_consts();

  void calculate_rot_matrix_miller(std::vector<int> vec_x_mil, std::vector<int> vec_y_mil);

  void set_xyz_mil_direction(std::string dir, int h, int k, int i, int l );


  static WzRotatedCrystal* create(const ModelOptions& options);

 protected:

  virtual void read_database ( );

  virtual void do_init(void) ;


  virtual void copy_from (const PhysicalModelInterface *rhs) ;


  virtual void do_init_alloy (const PhysicalModelInterface *comp_A, const PhysicalModelInterface *comp_B, double xa) ;


  virtual PhysicalModelInterface* create_new(void) const ;



 private:

  double a_lat;

  double c_lat;

};



inline WzRotatedCrystal* WzRotatedCrystal::create(const ModelOptions& options)
{
  return new WzRotatedCrystal(options);
}

#endif
