// $Id: ThermalSurfaceResistance.h 2069 2010-09-08 18:08:39Z gromano $

#ifndef _THERMALSURFACERESISTANCE_H_
#define _THERMALSURFACERESISTANCE_H_

#include "ThermalBoundaryModel.h"

#include "point.h"
#include "tensor_value.h"
#include "vector_value.h"
#include "tiber_dll.h"


namespace libMesh
{
  class Elem;
}



//! The base class for Poisson boundary conditions
class TBDLLOCAL ThermalSurfaceResistance : public ThermalBoundaryModel
{

  public:

    //! Destructor
    ~ThermalSurfaceResistance(void) {};

    //! Creator function
    static ThermalSurfaceResistance* create(const ModelOptions& options);

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


  private:

    //! Constructor
    ThermalSurfaceResistance(const ModelOptions& options);

  double _temperature;
  double _resistance;
};



inline
ThermalSurfaceResistance::ThermalSurfaceResistance(const ModelOptions& options) :
  ThermalBoundaryModel(options),
  _temperature(0),
  _resistance(0)
{
}



inline
ThermalSurfaceResistance*
ThermalSurfaceResistance::create(const ModelOptions& options)
{
  return new ThermalSurfaceResistance(options);
}





#endif // _POISSONDIRICHLET_H_
