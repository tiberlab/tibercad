// $Id$

#ifndef _PHYSICALMODEL_H_
#define _PHYSICALMODEL_H_

#include "PhysicalModelInterface.h"


//! Base class for basic physical models
/*!
 * A basic physical model is e.g. drift-diffusion, EFA,
 * \f$\mathbf k\cdot\mathbf p\f$ and similar.
 *
 * The models derived from this class will be usable as models
 * in the model list of the \c Material class.
 */
class PhysicalModel : public PhysicalModelInterface
{

  public:
    
    //! Destructor
    virtual ~PhysicalModel(void) {};


  protected:

    //! We don't want this to be instantiated directly
    PhysicalModel(void);
    
  private:
    
};

inline
PhysicalModel::PhysicalModel(void)
{
}



#endif // _PHYSICALMODEL_H_
