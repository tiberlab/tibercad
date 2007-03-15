#ifndef _THERMOELECTRICPOWER_H_
#define _THERMOELECTRICPOWER_H_

#include "PhysicalModelInterface.h"




//! Class to return the thermolectric power
/*!

The thermoelectric power must be in [mV/K]

*/

class ThermoelectricPower : public PhysicalModelInterface

{

public:
  
  //!Constructor 
  ThermoelectricPower()
  {  
   _eTEpower=0.0;
   _hTEpower=0.0;
  }

   //!Destructor
  ~ThermoelectricPower(){}


  //!provides thermoelectric power in simulation system [mV/K]
  inline void get_thermoelectric_power_e(double& TEpower) const; 

  inline void get_thermoelectric_power_h(double& TEpower) const; 
  
  
  inline  static ThermoelectricPower* create();
private:

 

protected:

  

  virtual void do_init (void);

  virtual void copy_from(const PhysicalModelInterface *rhs);

  virtual void read_database(void);

  virtual void read_bowing_parameters(void){};

  virtual void calculate_VCA (const PhysicalModelInterface *comp_A, const PhysicalModelInterface *comp_B, double xa); 


  inline  virtual PhysicalModelInterface* create_new (void) const;
 
 

 //!Thermoelectric power in simulation system. Units mV/(K)

  double _eTEpower;
  double _hTEpower;

 
 
};

void  ThermoelectricPower::get_thermoelectric_power_e(double& TEpower) const
{
  TEpower = _eTEpower;
}

void  ThermoelectricPower::get_thermoelectric_power_h(double& TEpower) const
{
  TEpower = _hTEpower;
}




ThermoelectricPower* ThermoelectricPower::create()
{
  return (new ThermoelectricPower());
}


PhysicalModelInterface* ThermoelectricPower::create_new (void) const
{
  return (new ThermoelectricPower() ); 
}
#endif
