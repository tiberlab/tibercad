// $Id$

#ifndef _LATTICEMISMATCH_H_
#define _LATTICEMISMATCH_H_

#include "BodyForceModel.h"

#include "point.h"
#include "tensor_value.h"
#include "vector_value.h"
#include "tiber_dll.h"

class Elem;



//! The base class for Poisson boundary conditions
class TBDLLOCAL LatticeMismatch : public BodyForceModel
{

  public:
 
  //! Destructor
  ~LatticeMismatch(void){};
  
  //! Creator function
  static LatticeMismatch* create(const ModelOptions& options);
  
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


    //! Create a new object of the same type
    virtual PhysicalModelInterface* create_new(void) const;


  private:
  

  //! Constructor
     LatticeMismatch(const ModelOptions& options);
  
};


inline
PhysicalModelInterface*
LatticeMismatch::create_new(void) const
{
  return new   LatticeMismatch(get_options());
}

inline
LatticeMismatch*
LatticeMismatch::create(const ModelOptions& options)
{ 
  return new  LatticeMismatch(options);
}




#endif // _GRAYMODEL_H_
