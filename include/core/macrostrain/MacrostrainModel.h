#ifndef _MACROSTRAINMODEL_H_
#define _MACROSTRAINMODEL_H_

#include "PhysicalModel.h"
#include "PhysicalModel.h"
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

  //!Destructor
  ~MacrostrainModel();


  //! creates a new object
  static MacrostrainModel* create(const ModelOptions& options);

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

  //! ID of the electric field solution
  ID id_E;

 protected:

  //!Constructor
  MacrostrainModel(const ModelOptions& options);

  virtual PhysicalModel* create_new (void) const;


  virtual void do_init();

  virtual void prepare_submodels();


};


inline MacrostrainModel* MacrostrainModel::create(const ModelOptions& options)
{
  return new  MacrostrainModel(options);
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
