// $Id: Clamp.h 2124 2010-10-22 14:00:17Z gromano $

#ifndef _BUILTINSTRAIN_H_
#define _BUILTINSTRAIN_H_

#include "BodyForceModel.h"

#include "tensor_value.h"




//! The base class for Poisson boundary conditions
class TBDLLOCAL BuiltInStrain : public BodyForceModel
{

  public:

    //! Destructor
    ~BuiltInStrain(void) {};

    //! Creator function
    static BuiltInStrain* create(const ModelOptions& options);


  protected:

    //! Initialize
  virtual void do_init(void);

    /* In some cases it might be useful to reimplement this: */
    // virtual void do_init_interface(const PhysicalModel* comp_A,
    //         const PhysicalModel* comp_B);


    /* This is not used here: */
    // virtual void read_database(void);


    /* We do not use this here: */
    // virtual void read_interface_database(void);



  private:

    //! Constructor
    BuiltInStrain(const ModelOptions& options);

   
};



inline
BuiltInStrain::BuiltInStrain(const ModelOptions& options) :
  BodyForceModel(options)
{
}



inline
BuiltInStrain*
BuiltInStrain::create(const ModelOptions& options)
{
  return new BuiltInStrain(options);
}




#endif // _POISSONDIRICHLET_H_
