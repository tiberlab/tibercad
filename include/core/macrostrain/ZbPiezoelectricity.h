#ifndef _ZBPIEZOELECTRICITY_H
#define _ZBPIEZOELECTRICITY_H

#include "Piezoelectricity.h"

//! Piezoelectric properties of a ZincBlende crystal 

class ZbPiezoelectricity : public Piezoelectricity
{
 public:


  //!Empty constructor
  ZbPiezoelectricity();

  //! constructor that sets module
  ZbPiezoelectricity(double e14);


  //! method that sets module
  void set_piezo_module(double e14);
  

  virtual void read_database (const Dummy &db);

  
  Tensor1  get_polariz_cryst(Tensor2Sym& strain_cryst);



 private:

  //!piezoelectric constant
  double e14;

};

#endif
