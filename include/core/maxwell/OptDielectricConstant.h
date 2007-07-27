#ifndef _OPTDIELECTRICCONSTANT_H_
#define _OPTDIELECTRICCONSTANT_H_

#include "PhysicalModelInterface.h"
#include "SimulationOptions.h"



//! Class to return the Optical Dielectric  constant
/*!

See lattice thermal conductivity.........

*/

class  OptDielectricConstant: public PhysicalModelInterface

{

 public:
  //!Constructor
  /*!

  A real and an imaginary tensor are  initialized.

  */
  OptDielectricConstant()
  {  
    _dielectric_constant_real=Tensor2Sym(0);
    _dielectric_constant_imag=Tensor2Sym(0);

  };

  //!Destructor
  ~OptDielectricConstant(){ } ;


  //!provides real part of dielectric constant   in simulation system
  void get_dielectric_real(Tensor2Sym& dielectric_constant_real) const; 

  //!provides imag. part of dielectric constant  in simulation system
  void get_dielectric_imag(Tensor2Sym& dielectric_constant_imag) const; 

 
     
  virtual void update_tensor(void)=0;


 private:


 protected:

  //! initialization
  virtual void do_init (void)=0;

  //! necessary to  assemble  alloys  materials
  virtual void copy_from(const PhysicalModelInterface *rhs);

  //! read material  data from  database
  virtual void read_database(void)=0;

  //!  read  bowing parameters of  material  present in  an  alloy
  virtual void read_bowing_parameters(void){};

  //! calculates dielectric  constant for  an  alloy, given the component materials  and their molar fraction.
  virtual void calculate_VCA (const PhysicalModelInterface *comp_A, const PhysicalModelInterface *comp_B, double xa); 

  //!Create a new model of the same type.
  virtual PhysicalModelInterface* create_new (void) const =0;


  //! real dielectric  tensor in simulation system. 
  Tensor2Sym _dielectric_constant_real;

  //! imag. dielectric  tensor in simulation system. 
  Tensor2Sym _dielectric_constant_imag;

  //!rotates  dielectric  tensor into the simulation system
  void rotate_to_calculation_system(const Tensor2Gen& RotMatrix);


};


inline
void
OptDielectricConstant::get_dielectric_real(Tensor2Sym& dielectric_constant_real) const
{

  dielectric_constant_real = _dielectric_constant_real;
  
}


inline
void
OptDielectricConstant::get_dielectric_imag(Tensor2Sym& dielectric_constant_imag) const
{

  dielectric_constant_imag = _dielectric_constant_imag;
  
}


inline
void 
OptDielectricConstant::rotate_to_calculation_system(const Tensor2Gen& RotMatrix)
{

  // generates dielectric  matrix in calculation system
  _dielectric_constant_real  = sym(RotMatrix * (_dielectric_constant_real * (RotMatrix.transpose())));
  _dielectric_constant_imag  = sym(RotMatrix * (_dielectric_constant_imag * (RotMatrix.transpose())));


}



#endif
