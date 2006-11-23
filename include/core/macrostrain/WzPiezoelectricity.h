#ifndef _WZPIEZOELECRICITY_H_
#define _WZPIEZOELECRICITY_H_

#include "Piezoelectricity.h"
//! Piezoelectric properties of wurtzite crystal
class WzPiezoelectricity : public Piezoelectricity
{
 public:

  WzPiezoelectricity();

  WzPiezoelectricity(double  e33, double e31, double e15, double Pz);

  void set_moduli(double  e33, double e31, double e15, double Pz);
  
  Tensor1 get_polariz_cryst(Tensor2Sym& strain_cryst);

  virtual void read_database (const Dummy &db);

 private:


  double e33;
  double e31;
  double e15;
  double Pz;
    
  

};

#endif
