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
  

   _eQfermi = 0.0;
   _hQfermi = 0.0;
   _Ec = 0.0;
   _Ev = 0.0;
   _e_mobility_term = 0.0;
   _h_mobility_term = 0.0;
   _TEmodel = "constant";
  }

   //!Destructor
  ~ThermoelectricPower(){}

  //! set the electro-chemical potential for electrons and holes
  void set_fermi_potential(double eQfermi, double hQfermi);

  //! set the band edge of conduction band and valence band
  void set_band_edges(double Ec, double Ev);

  //! set the local temperature
  void set_temperature(double Tloc);

   //! set the mobility term  for electron \f$\alpha_n\f$ and hole \f$\alpha_p\f$ 
  void set_mobility_term(double e_mobility_term, double h_mobility_term);


   //!provides electrons thermoelectric power [V/K]
  double  ThermoelectricPower::get_electrons_thermoelectric_power(void) const;

   //!provides holes thermoelectric power [V/K]
  double  ThermoelectricPower::get_holes_thermoelectric_power(void) const;
  
  //! Update the value of the thermoelectric power
  void re_init(void);
  
  static ThermoelectricPower* create();

private:

 
  double _Tloc; 

  double _eQfermi;
 
  double _hQfermi;

  double _Ec;

  double _Ev;

  double _e_mobility_term;  

  double _h_mobility_term;  

  //! Model for thermoelectric power 
  std::string _TEmodel;
 

protected:

  

  virtual void do_init (void);

  virtual void copy_from(const PhysicalModelInterface *rhs);

  virtual void read_database(void);

  virtual void read_bowing_parameters(void){};

  virtual void calculate_VCA (const PhysicalModelInterface *comp_A, const PhysicalModelInterface *comp_B, double xa); 


  virtual PhysicalModelInterface* create_new (void) const;
 
  //!Thermoelectric power in simulation system. Units V/K

  double _eTEpower;

  double _hTEpower;

 
 
};



inline
double  ThermoelectricPower::get_electrons_thermoelectric_power(void) const
{
  return _eTEpower;
}


inline
double  ThermoelectricPower::get_holes_thermoelectric_power(void) const
{
  return _hTEpower;
}


inline
void ThermoelectricPower::set_fermi_potential(double eQfermi, double hQfermi)
{
  _eQfermi = eQfermi;
  _hQfermi = hQfermi;
}


inline
void ThermoelectricPower::set_mobility_term(double e_mobility_term, double h_mobility_term)
{
  _e_mobility_term = e_mobility_term;

  _h_mobility_term = h_mobility_term;
 
}

inline
void ThermoelectricPower::set_band_edges(double Ec, double Ev)
{
  _Ec = Ec;
  _Ev = Ev;
}

inline
void ThermoelectricPower::set_temperature(double Tloc)
{
  _Tloc = Tloc;
}

inline
ThermoelectricPower* ThermoelectricPower::create()
{
  return (new ThermoelectricPower());
}

inline
PhysicalModelInterface* ThermoelectricPower::create_new () const
{
  return (new ThermoelectricPower() ); 
}




#endif
