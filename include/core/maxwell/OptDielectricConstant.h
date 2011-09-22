// $Id$

#ifndef _OPTDIELECTRICCONSTANT_H_
#define _OPTDIELECTRICCONSTANT_H_

#include "PhysicalModelInterface.h"
#include "SimulationOptions.h"



//! Class to return the Optical Dielectric  constant
/*!
  Used for Maxwell equations
*/

class  OptDielectricConstant: public PhysicalModelInterface

{

 public:
  //!Constructor
  /*!

  A real and an imaginary tensor are  initialized.

  */
  OptDielectricConstant(const ModelOptions& options) :
    PhysicalModelInterface(options)
  {  
    _dielectric_constant_real=Tensor2Sym(0);
    _dielectric_constant_imag=Tensor2Sym(0);

  };

  //!Destructor
  virtual ~OptDielectricConstant(){ } ;


  //!provides real part of dielectric constant   in simulation system
  void get_dielectric_real(Tensor2Sym& dielectric_constant_real) const; 

  //!provides imag. part of dielectric constant  in simulation system
  void get_dielectric_imag(Tensor2Sym& dielectric_constant_imag) const; 

  //!creates new model
  static OptDielectricConstant* create(const Material* mat, const ModelOptions &options);


 private:


 protected:

  //! initialization
  virtual void do_init (void)=0;

  //! read material  data from  database
  virtual void read_database(void)=0;

  //!  read  bowing parameters of  material  present in  an  alloy
  virtual void read_database_alloy(void){};

  //! calculates dielectric  constant for  an  alloy, given the component materials  and their molar fraction.
  virtual void do_init_alloy (const PhysicalModelInterface *comp_A, const PhysicalModelInterface *comp_B, double xa); 

  //!Create a new model of the same type.
  virtual PhysicalModelInterface* create_new (void) const =0;

  //! real dielectric  tensor in simulation system. 
  Tensor2Sym _dielectric_constant_real;

  //! imag. dielectric  tensor in simulation system. 
  Tensor2Sym _dielectric_constant_imag;

  //!rotates  dielectric  tensor into the simulation system
  void rotate_to_calculation_system(const Tensor2Gen& RotMatrix);

  //!name of the dielectric function model
  std::string _eps_model;


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

inline OptDielectricConstant* OptDielectricConstant::create(const Material* mat, const ModelOptions& options )
{
  std::string structure = mat->get_structure();
  return dynamic_cast<OptDielectricConstant*>(PhysicalModelInterface::create("opt_dielectric_constant_" + structure, mat, options));
}

#endif
