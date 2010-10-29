// $Id$

#ifndef _ANISOTROPICSTIFFNESS_H_
#define _ANISOTROPICSTIFFNESS_H_

#include "StiffnessModel.h"

#include "point.h"
#include "tensor_value.h"
#include "vector_value.h"
#include "tiber_dll.h"

class Elem;



//! The base class for Poisson boundary conditions
class TBDLLOCAL AnisotropicStiffness : public StiffnessModel
{

  public:
 
  //! Destructor
  ~AnisotropicStiffness(void) {};
  
  //! Creator function
  static AnisotropicStiffness* create(const ModelOptions& options);
  
  virtual void calculate(const Elem* elem, const Point& point){};

  protected:

    //! Initialize
    virtual void do_init(void);

    /* In some cases it might be useful to reimplement this: */
    // virtual void do_init_interface(const PhysicalModelInterface* comp_A,
    //         const PhysicalModelInterface* comp_B);


    /* This is not used here: */
  virtual void read_database(void);


    /* We do not use this here: */
    // virtual void read_interface_database(void);


    //! Create a new object of the same type
    virtual PhysicalModelInterface* create_new(void) const;


  private:
  
  double _c11;
  double _c12;
  double _c13;
  double _c33;
  double _c44;
  

  //! Constructor
    AnisotropicStiffness(const ModelOptions& options);
  
};


inline
PhysicalModelInterface*
AnisotropicStiffness::create_new(void) const
{
  return new   AnisotropicStiffness(get_options());
}

inline
AnisotropicStiffness*
AnisotropicStiffness::create(const ModelOptions& options)
{ 
  return new  AnisotropicStiffness(options);
}




#endif // _GRAYMODEL_H_
