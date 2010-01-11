// $Id$

#ifndef _PHYSICALMODEL_H_
#define _PHYSICALMODEL_H_

#include "PhysicalModelInterface.h"


//! Base class for basic physical models
/*!
 * A basic physical model is e.g. drift-diffusion, EFA,
 * \f$\mathbf k\cdot\mathbf p\f$ and similar.
 *
 * Often a model derived from this class will contain submodels.
 * They are created from the \c physical_model blocks in the input file.
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
    PhysicalModel(const ModelOptions& options);
    

  private:

    
};

inline
PhysicalModel::PhysicalModel(const ModelOptions& options)
 : PhysicalModelInterface(options)
{
}



#endif // _PHYSICALMODEL_H_
