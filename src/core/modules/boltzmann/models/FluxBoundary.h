// $Id: FluxBoundary.h 2069 2010-09-08 18:08:39Z gromano $

#ifndef _FLUXBOUNDARY_H_
#define _FLUXBOUNDARY_H_

#include "BoltzmannBoundaryModel.h"

#include "point.h"
#include "tensor_value.h"
#include "vector_value.h"
#include "tiber_dll.h"




//! The base class for Poisson boundary conditions
class TBDLLOCAL FluxBoundary : public BoltzmannBoundaryModel
{

  public:

    //! Destructor
    ~FluxBoundary(void) {};

    //! Creator function
    static FluxBoundary* create(const ModelOptions& options);

 //! Calculate for a point on the given side
    virtual void calculate(const Elem* elem, unsigned int side,
        const Point& point);
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
    FluxBoundary(const ModelOptions& options);

  double _heat_flux;
};



inline
FluxBoundary::FluxBoundary(const ModelOptions& options) :
  BoltzmannBoundaryModel(options),
  _heat_flux(0)
{
}



inline
FluxBoundary*
FluxBoundary::create(const ModelOptions& options)
{
  return new FluxBoundary(options);
}



inline
PhysicalModelInterface*
FluxBoundary::create_new(void) const
{
  return new FluxBoundary(get_options());
}

#endif 
