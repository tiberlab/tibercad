#ifndef _PARTICLETHERMALCONDUCTIVITY_H_
#define _PARTICLETHERMALCONDUCTIVITY_H_

#include "PhysicalModelInterface.h"
#include "tensor.h"




//! This class computes the particle thermal conductivity
/*!

The particles thermal conductivity are in [W/(K * cm)] units and are given by the Wiedemann-Franz Law:


\f$ k_n = TL\sigma_n \f$

\f$ k_p = TL\sigma_p   \f$

where L is the Lorenz number

\f$ L = \frac{4 k^2_b}{\pi q^2} = 2.45 10^-8 \frac{W \Omega}{K^2} \f$

*/


class ParticleThermalConductivity : public PhysicalModelInterface
{

public:
  
  //!Constructor 
   ParticleThermalConductivity()
  {  

    Tensor2Sym _kappa_e = Tensor2Sym(0);
   
    Tensor2Sym _kappa_h = Tensor2Sym(0);
    
     _Tloc = 0.0;
    
     _kappa_e_x = 0.0;
     
     _kappa_h_x = 0.0;

  }

   //!Destructor
  ~ParticleThermalConductivity(){}

  //! set the electron conducibility;
  void set_electrons_conducibility(double sigma_e);

  //! set the electron conducibility;
  void set_holes_conducibility(double  sigma_h);

  //! set the local temperature
  void set_temperature(double Tloc);

  //!provides electrons thermoelectric power [W / (cm * K)]
  void  get_electrons_thermal_conductivity(Tensor2Sym& kappa_e) const;

  //!provides holes thermoelectric power [W / (cm * K)]
  void  get_holes_thermal_conductivity(Tensor2Sym& kappa_h) const;
  
  //! Update the value of the thermoelectric power
  void re_init(void);
  
  static ParticleThermalConductivity* create();

private:

 
  Tensor2Sym _kappa_e;

  Tensor2Sym _kappa_h;

  double _Tloc; 

  double _sigma_e;
 
  double _sigma_h;

  double _kappa_e_x;

  double _kappa_h_x;

  std::string _particle_kappa_model; 

 
protected:

  

  virtual void do_init (void);

  virtual void copy_from(const PhysicalModelInterface *rhs);

  virtual void read_database(void);

  virtual void read_bowing_parameters(void){};

  virtual void calculate_VCA (const PhysicalModelInterface *comp_A, const PhysicalModelInterface *comp_B, double xa); 

  virtual PhysicalModelInterface* create_new (void) const;
 
 
 
};



inline
void ParticleThermalConductivity::get_electrons_thermal_conductivity(Tensor2Sym& kappa_e) const
{
  kappa_e = _kappa_e;
}

inline
void ParticleThermalConductivity::get_holes_thermal_conductivity(Tensor2Sym& kappa_h) const
{
  kappa_h = _kappa_h;
}




inline
void  ParticleThermalConductivity::set_electrons_conducibility(double sigma_e)
{
  _sigma_e = sigma_e;
 
}


inline
void  ParticleThermalConductivity::set_holes_conducibility(double sigma_h)
{
  _sigma_h = sigma_h;
 
}

inline
void  ParticleThermalConductivity::set_temperature(double Tloc)
{
  _Tloc = Tloc;
}


inline
ParticleThermalConductivity*  ParticleThermalConductivity::create()
{
  return (new  ParticleThermalConductivity());
}

inline
PhysicalModelInterface*  ParticleThermalConductivity::create_new () const
{
  return (new  ParticleThermalConductivity() ); 
}




#endif
