#ifndef _LATTICETHERMALCONDUCTIVITY_H_
#define _LATTICETHERMALCONDUCTIVITY_H_

#include "PhysicalModelInterface.h"
#include "SimulationOptions.h"



//! Class to return the lattice thermal Conductivity
/*!

The lattice thermal conductivity must be in W/(cm K)

*/

class LatticeThermalConductivity : public PhysicalModelInterface

{

public:
  
  //!Constructor 
  LatticeThermalConductivity()
  {  
   _conductivity=Tensor2Sym(0);
  };

   //!Destructor
  ~LatticeThermalConductivity(){};


  //!provides conductivity in simulation system W/(cm K)
  inline void get_conductivity(Tensor2Sym& conductivity) const; 

  //Temperature

  double temperature;
     
   virtual void update_tensor(void)=0;


private:

 

protected:

  virtual void do_init (void)=0;

  virtual void copy_from(const PhysicalModelInterface *rhs);

  virtual void read_database(void)=0;

  virtual void read_bowing_parameters(void){};

  virtual void calculate_VCA (const PhysicalModelInterface *comp_A, const PhysicalModelInterface *comp_B, double xa); 

  virtual PhysicalModelInterface* create_new (void) const =0;

  //!conductivity tensor in simulation system. Units W/(cm K)
  Tensor2Sym _conductivity;


  //!rotates conductivity into the simulation system
 inline void rotate_to_calculation_system(const Tensor2Gen& RotMatrix);

};

void  LatticeThermalConductivity::get_conductivity(Tensor2Sym& conductivity) const
{
   conductivity = _conductivity;
}

 void LatticeThermalConductivity::rotate_to_calculation_system(const Tensor2Gen& RotMatrix)
{

  // generates conductivity matrix in calculation system
  _conductivity = sym(RotMatrix * ( _conductivity * (RotMatrix.transpose())));


}



#endif
