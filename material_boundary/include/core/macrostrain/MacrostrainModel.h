#ifndef _MACROSTRAINMODEL_H_
#define _MACROSTRAINMODEL_H_

#include "PhysicalModel.h"
#include "PhysicalModelInterface.h"
#include "Stiffness.h"
#include "Piezoelectricity.h"
#include "RotatedCrystal.h"
#include "MacrostrainModelInterface.h"
#include "tensor.h"
#include "elem.h"

//!Class that contains all the objects, necessary for Macrostrain solver

class MacrostrainModel: public MacrostrainModelInterface
{
 public:
  //!Constructor
  MacrostrainModel();

  //!Destructor
  ~MacrostrainModel();


  //! creates a new object
  static MacrostrainModel* create();

  //! define an object that stores Young moduli
  void add_stiffness(Stiffness*  st);

  
  //! define the piezoelectricity
  void add_piezo(Piezoelectricity* pz);


  inline Stiffness* get_stiffness(void);


  inline Piezoelectricity* get_piezo(void);

  //!  calculates \f$ \sigma_{jk} = d_{i,jk}E_i \f$
  void get_converse_piezo_stress(Tensor2Sym& eps, const Elem* element);

 


 private:

  Stiffness* stiffness;


  Piezoelectricity* piezo;


  //! A pointer to the Poissons solver
  SimulationInterface* poisson;
  
  //!copy constructor should not be used
  MacrostrainModel (const MacrostrainModel &  t) {};

  //! set of ID's of the electric field components
  std::set< ID > Poisson_variables_ID;

  //!  ID's of the electric field components
  ID id_Ex; 
  ID id_Ey;
  ID id_Ez;

 protected:

  
  virtual PhysicalModelInterface* create_new (void) const;


  virtual void do_init_alloy (const PhysicalModelInterface *comp_A, const PhysicalModelInterface *comp_B, double xa);

  virtual void do_init();
  

};


inline MacrostrainModel* MacrostrainModel::create()
{
  return new  MacrostrainModel();
}


inline Stiffness* MacrostrainModel::get_stiffness(void)
{
  return( stiffness );
}

inline Piezoelectricity* MacrostrainModel::get_piezo(void)
{
  return( piezo );
}




#endif 
