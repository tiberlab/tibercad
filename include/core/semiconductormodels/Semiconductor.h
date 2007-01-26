#ifndef _SEMICONDUCTOR_H_
#define _SEMICONDUCTOR_H_

//!  A general crystal semiconductor class.
/*!
     The class can read parameters database and convert them into kp parameters 
*/

#include "tensor.h"
#include <vector>
#include "PhysicalModelInterface.h"

#include "KPparameters.h"

class Semiconductor : public PhysicalModelInterface
{
 public:

  //!Constructor
   Semiconductor(void){};
 

  //!Desctructor
   virtual ~Semiconductor(void) {};




  //! Calculates k.p parameters in atomic units for 6 band valence band calculation
  virtual KPparams   calculate_6x6_kp_params (void )=0;

  //! Calculates k.p parameters in atomic units for 8 band valence band calculation
  virtual KPparams   calculate_8x8_kp_params (void )=0;


  //! Calculates k.p parameters in atomic units for 6 or 8 band valence band calculation
  KPparams   calculate_kp_params (std::string kp_model );


  
  //! creates new object
  static Semiconductor* create(const std::string& name,  const ModelOptions& options);
 

 private:

  //! Hartree energy in eV
  static const double Hartree;
   

 

 protected:

  virtual PhysicalModelInterface* create_new(void) const;

  virtual void copy_from (const PhysicalModelInterface *rhs) ;

  virtual void do_init (void);

  virtual void read_database(void) = 0;

  virtual void read_bowing_parameters(void) = 0;
 
  virtual void calculate_VCA (const PhysicalModelInterface *comp_A, const PhysicalModelInterface *comp_B, double xa);

};

inline Semiconductor* Semiconductor::create(const std::string& name,  const ModelOptions& options)
{
  return dynamic_cast<Semiconductor*> (PhysicalModelInterface::create("semicond_" + name ,options));
}
 


#endif
