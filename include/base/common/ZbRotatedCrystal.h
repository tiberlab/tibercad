#ifndef _ZBROTATEDCRYSTAL_H_
#define _ZBROTATEDCRYSTAL_H_

#include "RotatedCrystal.h"

class TBDLLOCAL ZbRotatedCrystal : public RotatedCrystal
{

 public:

  void set_lat_const(const double a);

  virtual void calculate_lat_consts();

  void calculate_rot_matrix_miller(std::vector<int> vec_x_mil, std::vector<int> vec_y_mil);

  void set_xyz_mil_direction(std::string dir, int h, int k, int l);


  static ZbRotatedCrystal* create(const ModelOptions& options);

 protected:

  ZbRotatedCrystal(const ModelOptions& options);

  virtual void read_database ( ) ;

 
  virtual void do_init(void) ;


  virtual void copy_from (const PhysicalModelInterface *rhs) ;


  virtual void do_init_alloy (const PhysicalModelInterface *comp_A, const PhysicalModelInterface *comp_B, double xa) ;


  virtual PhysicalModelInterface* create_new(void) const ;


  

 private:

  //! Lattice constant
  double a_lat; 


};

inline ZbRotatedCrystal* ZbRotatedCrystal::create(const ModelOptions& options)
{
  return new ZbRotatedCrystal(options);
}

#endif
