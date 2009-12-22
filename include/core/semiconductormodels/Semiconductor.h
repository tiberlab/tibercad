#ifndef _SEMICONDUCTOR_H_
#define _SEMICONDUCTOR_H_

//!  A general crystal semiconductor class.
/*!
     The class can read parameters database and convert them into kp parameters 
*/

#include "tensor.h"
#include <vector>
#include "PhysicalModelInterface.h"
#include "TemperatureInterface.h"

#include "KPparameters.h"

class Semiconductor : public PhysicalModelInterface
{
 public:

  //!Constructor
   Semiconductor(void);
 

  //!Desctructor
   virtual ~Semiconductor(void) {};




  //! Calculates k.p parameters in atomic units for 6 band valence band calculation
  virtual KPparams   calculate_6x6_kp_params (void );

  //! Calculates k.p parameters in atomic units for 8 band valence band calculation
  virtual KPparams   calculate_8x8_kp_params (void );



  //! Calculates k.p parameters in atomic units for 6 band valence band calculation with temperature dependence
  virtual KPparams   calculate_6x6_kp_params (const Elem* element, const Point& point) ;

  //! Calculates k.p parameters in atomic units for 8 band valence band calculation with temperature dependence
  virtual KPparams   calculate_8x8_kp_params (const Elem* element, const Point& point) ;


  //! apply varshni formulas
  virtual void apply_temperature(void) = 0;



  //! Calculates k.p parameters in atomic units for 6 or 8 band valence band calculation
  KPparams   calculate_kp_params (std::string kp_model );


  
  //! creates new object
  static Semiconductor* create(const std::string& name,  const ModelOptions& options);
 
  //! sets temperature
  void set_temperature(double T);


  //! sets temperature from a simulation
  void set_temperature(const Elem* element, const Point& point);


 private:

  //! Hartree energy in eV
  static const double Hartree;
   
  
 protected:

  virtual PhysicalModelInterface* create_new(void) const = 0;

  virtual void do_init (void);

  virtual void read_database(void) = 0;

  virtual void read_database_alloy(void) = 0;
 
  virtual void do_init_alloy (const PhysicalModelInterface *comp_A, const PhysicalModelInterface *comp_B, double xa);


  virtual void do_do_init_alloy (const PhysicalModelInterface *comp_A, const PhysicalModelInterface *comp_B, double xa) = 0;


  //! if true band gap is temperature 
  bool _consider_temperature;
  
  
  //! lattice temperature [K]
  double _temperature;

  //! interface for temperature acquisition
  TemperatureInterface temp_interface;


  //! Calculates k.p parameters in atomic units for 6 band valence band calculation
  virtual KPparams   do_calculate_6x6_kp_params (void )=0;

  //! Calculates k.p parameters in atomic units for 8 band valence band calculation
  virtual KPparams   do_calculate_8x8_kp_params (void )=0;

  //! pointer to a parent model A (needed only for alloys) 
  const Semiconductor* modelA;

  //! pointer to a parent model B (needed only for alloys) 
  const Semiconductor* modelB;

  //! molar fraction (a local copy)
  double _xa;


};

inline 
Semiconductor* Semiconductor::create(const std::string& name,  const ModelOptions& options)
{
  return dynamic_cast<Semiconductor*> (PhysicalModelInterface::create("semicond_" + name, options));
}
 


inline 
void  Semiconductor::set_temperature(double T)
{
  _temperature = T;
}


inline 
void Semiconductor::set_temperature(const Elem* element, const Point& point )
{
   _temperature = temp_interface.get_temperature (element, point);
}

inline  
KPparams   Semiconductor::calculate_6x6_kp_params()
{
  if (_consider_temperature) apply_temperature();

  return(do_calculate_6x6_kp_params());
}

inline 
KPparams   Semiconductor::calculate_8x8_kp_params()
{ 

  if (_consider_temperature) apply_temperature();

  return(do_calculate_8x8_kp_params());
}




inline
KPparams   Semiconductor::calculate_6x6_kp_params (const Elem* element, const Point& point) 
{

  _temperature = temp_interface.get_temperature ( element, point);
  return(calculate_6x6_kp_params());

}

  
inline  
KPparams  Semiconductor:: calculate_8x8_kp_params (const Elem* element, const Point& point) 
{
  _temperature = temp_interface.get_temperature ( element, point);

  return(calculate_8x8_kp_params());

}




#endif
