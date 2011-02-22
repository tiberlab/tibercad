// $Id$

#ifndef _ISOTROPICSTIFFNESS_H_
#define _ISOTROPICSTIFFNESS_H_

#include "StiffnessModel.h"

#include "point.h"
#include "tensor_value.h"
#include "vector_value.h"
#include "tiber_dll.h"

class Elem;



//! The base class for Poisson boundary conditions
class TBDLLOCAL IsotropicStiffness : public StiffnessModel
{

  public:
 
  //! Destructor
  ~IsotropicStiffness(void) {};
  
  //! Creator function
  static IsotropicStiffness* create(const ModelOptions& options);
  
  virtual void calculate(const Elem* elem, const Point& point){};

  protected:

    //! Initialize
    virtual void do_init(void);

    /* In some cases it might be useful to reimplement this: */
    // virtual void do_init_interface(const PhysicalModelInterface* comp_A,
    //         const PhysicalModelInterface* comp_B);


    /* This is not used here: */
  // virtual void read_database(void);


    /* We do not use this here: */
    // virtual void read_interface_database(void);



  private:
  
  double _young;
  double _poisson;
  

  //! Constructor
    IsotropicStiffness(const ModelOptions& options);
  
};




inline
IsotropicStiffness*
IsotropicStiffness::create(const ModelOptions& options)
{ 
  return new  IsotropicStiffness(options);
}




#endif // _GRAYMODEL_H_
