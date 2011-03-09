// $Id: HeatReservoir.h 2075 2010-09-15 09:40:22Z gromano $

#ifndef _HEATRESERVOIR_H_
#define _HEATRESERVOIR_H_

#include "ThermalBoundaryModel.h"

#include "point.h"
#include "tensor_value.h"
#include "vector_value.h"
#include "tiber_dll.h"


class Elem;


//! The base class for Poisson boundary conditions
class TBDLLOCAL HeatReservoir : public ThermalBoundaryModel
{

  public:

    //! Destructor
    ~HeatReservoir(void) {};

    //! Creator function
    static HeatReservoir* create(const ModelOptions& options);

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
    HeatReservoir(const ModelOptions& options);

  double _temperature;
};



inline
HeatReservoir::HeatReservoir(const ModelOptions& options) :
  ThermalBoundaryModel(options),
  _temperature(0)
{
}



inline
HeatReservoir*
HeatReservoir::create(const ModelOptions& options)
{
  return new HeatReservoir(options);
}



inline
PhysicalModelInterface*
HeatReservoir::create_new(void) const
{
  return new HeatReservoir(get_options());
}

#endif // _POISSONDIRICHLET_H_
