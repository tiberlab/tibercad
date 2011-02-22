// $Id$

#ifndef _CONSTANTBODYFORCE_H_
#define _CONSTANTBODYFORCE_H_

#include "BodyForceModel.h"

#include "point.h"
#include "tensor_value.h"
#include "vector_value.h"
#include "tiber_dll.h"

class Elem;



//! The base class for Poisson boundary conditions
class TBDLLOCAL ConstantBodyForce : public BodyForceModel
{

  public:
 
  //! Destructor
  ~ConstantBodyForce(void) {};
  
  //! Creator function
  static ConstantBodyForce* create(const ModelOptions& options);
  
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
  

  //! Constructor
     ConstantBodyForce(const ModelOptions& options);
  
};




inline
ConstantBodyForce*
ConstantBodyForce::create(const ModelOptions& options)
{ 
  return new  ConstantBodyForce(options);
}




#endif // _GRAYMODEL_H_
