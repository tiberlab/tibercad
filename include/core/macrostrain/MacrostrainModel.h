#ifndef _MACROSTRAINMODEL_H_
#define _MACROSTRAINMODEL_H_

#include "PhysicalModel.h"
#include "PhysicalModelInterface.h"
#include "Stiffness.h"
#include "Piezoelectricity.h"
#include "RotatedCrystal.h"
#include "MacrostrainModelInterface.h"

//!Class that contains all the objects, necessary for Macrostrain solver

class MacrostrainModel: public MacrostrainModelInterface
{
 public:
  //!Constructor
  MacrostrainModel();

  //!Destructor
  ~MacrostrainModel();

  static MacrostrainModel* create();

  //! define an object that stores Young moduli
  void add_stiffness(Stiffness*  st);

  
  //! define the piezoelectricity
  void add_piezo(Piezoelectricity* pz);


  inline Stiffness* get_stiffness(void);


  inline Piezoelectricity* get_piezo(void);


 


 private:

  Stiffness* stiffness;


  Piezoelectricity* piezo;


 
  
  //!copy constructor should not be used
  MacrostrainModel (const MacrostrainModel &  t) {};
 

 protected:

  
  virtual PhysicalModelInterface* create_new (void) const;


  virtual void copy_from(const PhysicalModelInterface *rhs);


  virtual void read_database (void);
 

  virtual void read_bowing_parameters (void) {};


  virtual void calculate_VCA (const PhysicalModelInterface *comp_A, const PhysicalModelInterface *comp_B, double xa);

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
