#ifndef _SEMICONDUCTOR_H_
#define _SEMICONDUCTOR_H_


#include "PhysicalModelInterface.h"
#include "libMeshDefs.h"
#include "TemperatureInterface.h"
#include "KPparameters.h"

//#include "tensor.h"
//#include "elem.h"
#include <vector>

//! A general crystal semiconductor class.
/*!
 *  The class can read parameters database and convert them into kp parameters
 */
class TBDLEXPORT Semiconductor : public PhysicalModelInterface
{
 public:
 

  //!Desctructor
   virtual ~Semiconductor(void) {};


  //! apply varshni formulas
  virtual void apply_temperature(void) = 0;



  //! Calculates k.p parameters in atomic units for 6 or 8 band valence band calculation
  void calculate_kp_params(KPparams& par);

  //  ! used to set or reset kp model
  //void set_kp_model(std::string kp_model);
  
  //! creates new object
  static Semiconductor* create(const Material* mat,  const ModelOptions& options);
 
  //! sets temperature
  void set_temperature(double T);


  //! sets temperature from a simulation
  void set_temperature(const Elem* element, const Point& point);

  
 protected:

  //!Constructor
   Semiconductor(const ModelOptions& options);

  virtual PhysicalModelInterface* create_new(void) const = 0;

  virtual void do_init (void);

  virtual void read_database(void) = 0;

  virtual void do_init_alloy (const PhysicalModelInterface *comp_A, const PhysicalModelInterface *comp_B, double xa);

  //! if true band gap is temperature 
  bool _consider_temperature;
  
  bool _couple_bands;
  
  
  //! lattice temperature [K]
  double _temperature;

  //! interface for temperature acquisition
  TemperatureInterface temp_interface;

  //! Calculates k.p parameters in atomic units (will depend on _kp_model)
  virtual void do_calculate_kp_params (KPparams& par)=0;

  //! pointer to a parent model A (needed only for alloys) 
  const Semiconductor* modelA;

  //! pointer to a parent model B (needed only for alloys) 
  const Semiconductor* modelB;

  //! molar fraction (a local copy)
  double _xa;

  std::string _kp_model;

  std::string _spurious;

};




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

//inline
//void Semiconductor::set_kp_model(std::string kp_model)
//{
//  _kp_model = kp_model;
//}


inline  
void Semiconductor::calculate_kp_params(KPparams& par)
{
  if (_consider_temperature) apply_temperature();

  do_calculate_kp_params(par);
}





#endif
