#ifndef _POISSONMODEL_H_
#define _POISSONMODEL_H_


#include "PhysicalModel.h"
#include "elem.h"
#include "SimulationInterface.h"
#include "ChargeDensityModel.h"
#include "DielectricModel.h"

//!Class that contains all the physical quantities necessary for the POISSON solver
class PoissonModel: public PhysicalModel
{
 public:

  //!Constructor
  PoissonModel(void);
  
  //!Destructor
  ~PoissonModel(void);
  
  //! creates a new object
  static  PoissonModel* create(void);
  
  void 	re_init(void);
  
  //! update the charge density
  void  update_charge_density(void);
  
  //!Set the current element
  void set_element(const Elem* elem);
  
  //!Return the charge density for the current element
  double get_charge_density(); 
  
  //!Return the charge density for the current element
  Tensor2Sym get_dielectric_constant(); 



 private:
    
   //! Current element 
   const Elem* _elem; 

   //!copy constructor should not be used
   PoissonModel (const PoissonModel&  t) {};

   //!relative dielectric constant
   Tensor2Sym _epsilon;

   //!charge density (electron/cm^3)
   double _charge_density;

   //!A pointer to cherge density object
   ChargeDensityModel* chd_model;

   //!A pointer to dielectric constant model
   DielectricModel* dielectric_model;


  
 protected:

  virtual PhysicalModelInterface* create_new (void) const; 


  virtual void copy_from(const PhysicalModelInterface *rhs){};


  virtual void read_database (void){};
 

  virtual void read_bowing_parameters (void) {};


  virtual void calculate_VCA (const PhysicalModelInterface *comp_A, const PhysicalModelInterface *comp_B, double xa);

  virtual void do_init();

};


inline 
PoissonModel* 
PoissonModel::create(void)
{
  return new PoissonModel;
}

inline
void 
PoissonModel::set_element(const Elem* elem)
{

  _elem = elem;

}


inline   
double
PoissonModel::get_charge_density()
{

 return _charge_density;

}



inline   
Tensor2Sym
PoissonModel::get_dielectric_constant()
{

 
  return _epsilon;

}

#endif
