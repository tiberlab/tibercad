#ifndef _ROTATED_CRYSTAL_H_
#define _ROTATED_CRYSTAL_H_

#include "tensor.h"
#include <cmath>
#include <vector>

#include "PhysicalModelInterface.h"
#include "PhysicalModel.h"

class RotatedCrystal : public PhysicalModelInterface
{

 public:

   //! Create a RotatedCrystal object
   static RotatedCrystal* create(const std::string& name,
       const ModelOptions& options);
   
  //!Rotation matrix  \f$ {\bf x}^{calc} = R {\bf x}^{cryst} \f$  \f$ 


  /*!In matrix notation:  
    vector transformation:    vec_calc = R * vec_crystal
    2D tensor transformation: Matr_calc = R * Mat_crystal * (R^T)
  */
  Tensor2Gen RotMatrix;

  //! calculates Bravais vectors, rotation matrix and lattice constants in calculation system 
  virtual void calculate_lat_consts()=0;


  //! Calculates rotation matrix
  /*!
    \param  vec_x direction in crystal system that is x direction in the calculation system 
    \param  vec_y direction in crystal system that is y direction in the calculation system 
  */
  void calculate_rot_matrix(const Tensor1& vec_x,const Tensor1& vec_y);



  //! Calculates rotation matrix basing on miller indexes 
  virtual void calculate_rot_matrix_miller(std::vector<int> vec_x_mil, std::vector<int> vec_y_mil) = 0;



  //! returns lattice constants along growth directions
  /*!
    \param lat_cont  at the exit contains lattice constants
   */
  void get_lat_const(double lat_cont[3]); 


  //! calculates lattice matching deformation tensor in calculation system
  /*!
    now we assume that the calculation system basis is parallel to conventinal lattice cell basis 
    \f$ \varepsilon_{ij}^{0} = \delta_{ij} \frac {a^{substrate}_{i} - a_{i}} {a_i}  \f$
    \param lat_cont_substrate substrate lattice constants (conventional cell)
   */
  Tensor2Sym get_eps0(double lat_cont_substrate[3]);


  //! calculates constant part of the lattice matching tensor \f$ \varepsilon^{0}_{ij} \f$
  /*!
    
    \param lat_cont_substrate substrate lattice constants (conventional cell)
    \param eps0_var_log this tensor defines which lattice constants are fixed.

  */
  Tensor2Sym get_const_eps0(double lat_cont_substrate[3], Tensor2Sym& eps0_var_log);

  //! calculates variable part of the lattice matching tensor
  /*!
    \param name of the additional variable
  */

  Tensor2Sym get_var_eps0(std::string  var_name);


 


 protected:

  //! Constructor
  RotatedCrystal();

  

  //! Miller indexes of the growth direction x. At the moment it has to be parallel to the x direction of the simulation system
  std::vector<int> x_miller; 
  

  //! Miller indexes of the growth direction y. At the moment it has to be parallel to the y direction of the simulation system
  std::vector<int> y_miller;


  //! Miller indexes of the growth direction z. At the moment it has to be parallel to the z direction of the simulation system
  std::vector<int> z_miller;


  //! lattice constants along the 3 main growth directions
  double  lat_const_calc[3];


  //!read lattice constant from the database
  virtual void read_database ( ) =0;


 
  virtual void do_init(void) = 0;


  virtual void copy_from (const PhysicalModelInterface *rhs) = 0;


  virtual void calculate_VCA (const PhysicalModelInterface *comp_A, const PhysicalModelInterface *comp_B, double xa) = 0;


  virtual PhysicalModelInterface* create_new(void) const = 0;
 


 private:


 



};

//
// inline members
//

inline
RotatedCrystal*
RotatedCrystal::create(const std::string& name, const ModelOptions& options)
{
  return dynamic_cast<RotatedCrystal*>(
      PhysicalModelInterface::create("cryst_" + name, options));
}


#endif
