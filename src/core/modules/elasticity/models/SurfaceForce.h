// $Id$

#ifndef _SURFACEFORCE_H_
#define _SURFACEFORCE_H_

#include "BoundaryModel.h"

#include "point.h"
#include "tensor_value.h"
#include "vector_value.h"
#include "TensorOperators.h"
#include "tiber_dll.h"

class Elem;


//! The base class for Poisson boundary conditions
class TBDLLOCAL SurfaceForce : public BoundaryModel
{

  public:

    //! Destructor
    ~SurfaceForce(void) {};

    //! Creator function
    static SurfaceForce* create(const ModelOptions& options);


    //! Calculate for a point on the given side
    virtual void calculate(const Elem* elem, unsigned int side,
			   const Point& point){};


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
    SurfaceForce(const ModelOptions& options);

    //! The boundary potential
    //double _potential;

};



inline
SurfaceForce::SurfaceForce(const ModelOptions& options) :
  BoundaryModel(options)
{
}



inline
SurfaceForce*
SurfaceForce::create(const ModelOptions& options)
{
  return new SurfaceForce(options);
}



inline
PhysicalModelInterface*
SurfaceForce::create_new(void) const
{
  return new SurfaceForce(get_options());
}

#endif // _POISSONDIRICHLET_H_
