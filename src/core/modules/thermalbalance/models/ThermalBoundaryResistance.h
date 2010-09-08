// $Id$

#ifndef _THERMALBOUNDARYRESISTANCE_H_
#define _THERMALBOUNDARYRESISTANCE_H_

#include "ThermalBoundaryModel.h"

#include "point.h"
#include "tensor_value.h"
#include "vector_value.h"
#include "tiber_dll.h"


class Elem;


//! The base class for Poisson boundary conditions
class TBDLLOCAL ThermalBoundaryResistance : public ThermalBoundaryModel
{

  public:

    //! Destructor
    ~ThermalBoundaryResistance(void) {};

    //! Creator function
    static ThermalBoundaryResistance* create(const ModelOptions& options);

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
    ThermalBoundaryResistance(const ModelOptions& options);

  double _temperature;
  double _resistance;
};



inline
ThermalBoundaryResistance::ThermalBoundaryResistance(const ModelOptions& options) :
  ThermalBoundaryModel(options),
  _temperature(0),
  _resistance(0)
{
}



inline
ThermalBoundaryResistance*
ThermalBoundaryResistance::create(const ModelOptions& options)
{
  return new ThermalBoundaryResistance(options);
}



inline
PhysicalModelInterface*
ThermalBoundaryResistance::create_new(void) const
{
  return new ThermalBoundaryResistance(get_options());
}

#endif // _POISSONDIRICHLET_H_
