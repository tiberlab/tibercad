// $Id$

#ifndef _WZPIEZOELECRICITY_H_
#define _WZPIEZOELECRICITY_H_

#include "Piezoelectricity.h"

//! Piezoelectric properties of wurtzite crystal
class WzPiezoelectricity : public Piezoelectricity
{
 public:


  void set_moduli(double  e33, double e31, double e15);
  
  //!calculates piezopolarization in crystal system
  /*!
    \param strain_cryst strain tensor in crystal system

    \f$ 
    P_x = 2e_{15} \epsilon_{xz} \\ 
    P_y = 2e_{15} \epsilon_{yz} \\
    P_z = e_{31}  \epsilon_{xx} +  e_{31}  \epsilon_{yy} + e_{33}  \epsilon_{zz}
    \f$
  */
  Tensor1 get_polariz_cryst(Tensor2Sym& strain_cryst);

  virtual void calculate_product_by_vector(const Tensor1& f, Tensor2Sym& r) const;


  static WzPiezoelectricity* create(const ModelOptions& options);

 private:


  double e33;
  double e31;
  double e15;


 protected:

  WzPiezoelectricity(const ModelOptions& options);

  virtual void read_database(void);

  virtual void do_init(void);


  virtual PhysicalModelInterface* create_new(void) const;
 

};


inline WzPiezoelectricity* WzPiezoelectricity::create(const ModelOptions& options)
{
  return new WzPiezoelectricity(options);
}

#endif
