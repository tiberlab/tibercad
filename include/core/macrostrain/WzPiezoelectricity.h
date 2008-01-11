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
