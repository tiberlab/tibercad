#ifndef _THERMOELECTRICPOWER_H_
#define _THERMOELECTRICPOWER_H_

#include "PhysicalModelInterface.h"




//! This class computes the thermolectric power
/*!

The thermoelectric power must be in [V/K]

  
 \f$ P_n = -\frac{k_b}{q}\left(\frac{5}{2} + \alpha_n + E_c + q \varphi_n   \right)  \f$

  \f$ P_p = +\frac{k_b}{q}\left(\frac{5}{2} + \alpha_p  - q \varphi_n - E_v \right)  \f$

  where
 
 \f$  \alpha = \frac{T}{\mu_n} \frac {\partial \mu_n}{\partial T}\f$
 \f$  \alpha = \frac{T}{\mu_p} \frac {\partial \mu_p}{\partial T}\f$

*/


class ThermoelectricPower : public PhysicalModelInterface

{

public:
  
  //!Constructor 
  ThermoelectricPower()
  {  
  _eTEpower=0.0;
  _hTEpower=0.0;
  KfracQ = 2.2944e-7;  // V/K

   _eQfermi = 0.0;
   _hQfermi = 0.0;
   _Ec = 0;
   _Ev = 0;
   _e_mobility_term = 0;
   _h_mobility_term = 0;
   _eTEmodel = "constant";

   _hTEmodel = "constant";
  }

   //!Destructor
  ~ThermoelectricPower(){}

  //! set the electro-chemical potential for electrons and holes
  inline void set_Qfermi(double& eQfermi, double& hQfermi);

  //! set the band edge of conduction band and valence band
  inline void set_band_edge(double& Ec, double& Ev);

   //! set the mobility term  for electron \f$\alpha_n\f$ and hole \f$\alpha_p\f$ 
  inline void set_mobility_term(double& e_mobility_term, double& h_mobility_term);

   //!provides electron thermoelectric power in simulation system [V/K]
  inline void get_thermoelectric_power_e(double& TEpower) const; 

  //!provides hole thermoelectric power in simulation system [V/K]
  inline void get_thermoelectric_power_h(double& TEpower) const; 
  
  //! Update the value of the thermoelectric power
  void update_tensor(void);
  
  inline  static ThermoelectricPower* create();

private:

  double KfracQ; 

  double _eQfermi;
 
  double _hQfermi;

  double _Ec;

  double _Ev;

  double _e_mobility_term;  

  double _h_mobility_term;  

  //! Model for electron thermoelectric power 
  std::string _eTEmodel;

  //! Model for hole thermoelectric power 
  std::string _hTEmodel;
 

protected:

  

  virtual void do_init (void);

  virtual void copy_from(const PhysicalModelInterface *rhs);

  virtual void read_database(void);

  virtual void read_bowing_parameters(void){};

  virtual void calculate_VCA (const PhysicalModelInterface *comp_A, const PhysicalModelInterface *comp_B, double xa); 


  inline  virtual PhysicalModelInterface* create_new (void) const;
 
  //!Thermoelectric power in simulation system. Units V/K

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

inline
void ThermoelectricPower::set_Qfermi(double& eQfermi, double& hQfermi)
{
  _eQfermi = eQfermi;
  _hQfermi = hQfermi;
}


inline
void ThermoelectricPower::set_mobility_term(double& e_mobility_term, double& h_mobility_term)
{
  _e_mobility_term = e_mobility_term;

  _h_mobility_term = h_mobility_term;
 
}

inline
void ThermoelectricPower::set_band_edge(double& Ec, double& Ev)
{
  _Ec = Ec;
  _Ev = Ev;
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
