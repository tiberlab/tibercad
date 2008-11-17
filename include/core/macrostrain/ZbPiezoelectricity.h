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
  

 

  //! calculate piezoelectric polarization in crystal system  
  Tensor1  get_polariz_cryst(Tensor2Sym& strain_cryst);


  virtual void calculate_product_by_vector(const Tensor1& f, Tensor2Sym& r) const;


 


  static ZbPiezoelectricity* create();


 private:

  //!piezoelectric constant
  double e14;

  //!bowing parameter for the piezoelectric constant
  double e14_bow;

 protected:

  virtual void read_database(void) ;

  virtual void read_database_alloy(void);


  virtual void do_init(void);


  virtual void do_init_alloy (const PhysicalModelInterface *comp_A, const PhysicalModelInterface *comp_B, double xa);

  
  virtual PhysicalModelInterface* create_new(void) const;



};


inline ZbPiezoelectricity*  ZbPiezoelectricity::create()
{
  return new ZbPiezoelectricity();
}

#endif
