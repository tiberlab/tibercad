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

  virtual void calculate_product_by_vector(const Tensor1& f, Tensor2Sym& r) const;


  static WzPiezoelectricity* create();

 private:


  double e33;
  double e31;
  double e15;
  double Pz;


  double e33_bow;
  double e31_bow;
  double e15_bow;
  double Pz_bow;
    

 protected:

  virtual void read_database ( ) ;


  virtual void do_init(void);


  virtual void copy_from (const PhysicalModelInterface *rhs);


  virtual void calculate_VCA (const PhysicalModelInterface *comp_A, const PhysicalModelInterface *comp_B, double xa);

  
  virtual PhysicalModelInterface* create_new(void) const;
 

};


inline WzPiezoelectricity* WzPiezoelectricity::create()
{
  return new WzPiezoelectricity();
}

#endif
