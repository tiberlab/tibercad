#ifndef _PIEZOELECTRICITY_H_
#define _PIEZOELECTRICITY_H_

#include "tensor.h"
#include "xtensor.h"
#include <cmath>
#include "PhysicalProperties.h"

class Piezoelectricity : public PhysicalProperties
{
 public:

  //!Empty constructor
  Piezoelectricity();

  //! returns polarization (piezo + pyro) in crystal system
  virtual Tensor1  get_polariz_cryst(Tensor2Sym& strain_cryst) = 0;

  
  virtual void read_database (const Dummy &db) {};

};


#endif
