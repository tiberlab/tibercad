#ifndef _HEATMODEL_H_
#define _HEATMODEL_H_


#include "PhysicalModel.h"
#include "LatticeThermalConductivity.h"
#include "ThermoelectricPower.h"



//!Class that contains all the object, necessary for Heat Transport solver
class HeatModel: public PhysicalModel
{
 public:
  //!Constructor
  HeatModel();

  //!Destructor
  ~HeatModel();
   
   //! creates a new object
  static HeatModel* create();
  
  inline const LatticeThermalConductivity* get_lattice_conductivity(void) const;

  inline const ThermoelectricPower* get_thermoelectric_power(void) const;

  
 private:

   LatticeThermalConductivity* kappa;

   ThermoelectricPower* thermoelectric_power;
  

  //!copy constructor should not be used
  HeatModel (const HeatModel &  t) {};
  
 protected:

  virtual PhysicalModelInterface* create_new (void) const;


  virtual void copy_from(const PhysicalModelInterface *rhs){};


  virtual void read_database (void){};
 

  virtual void read_bowing_parameters (void) {};


  virtual void calculate_VCA (const PhysicalModelInterface *comp_A, const PhysicalModelInterface *comp_B, double xa);

  virtual void do_init();

};



inline HeatModel* HeatModel::create()
{
  return new  HeatModel();
}


inline const LatticeThermalConductivity* HeatModel::get_lattice_conductivity(void) const
{
  return( kappa );
}

inline const ThermoelectricPower* HeatModel::get_thermoelectric_power(void) const
{
  return( thermoelectric_power);
}



#endif
