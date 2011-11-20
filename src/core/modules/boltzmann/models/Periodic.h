// $Id: HeatReservoir.h 2075 2010-09-15 09:40:22Z gromano $

#ifndef _PERIODC_H_
#define _PERIODC_H_

#include "BoltzmannBoundaryModel.h"

#include "point.h"
#include "tensor_value.h"
#include "vector_value.h"
#include "tiber_dll.h"


class Elem;

//! The base class for Poisson boundary conditions
class TBDLLOCAL Periodic : public BoltzmannBoundaryModel
{

  public:

    //! Destructor
    ~Periodic(void) {};

    //! Creator function
    static Periodic* create(const ModelOptions& options);

   //! Calculate for a point on the given side
    virtual void calculate(const Elem* elem, unsigned int side,
        const Point& point);

    //! Calculate for a point on the given side
    //virtual void get_periodicity(const Point& point);

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
    Periodic(const ModelOptions& options);

  double _deltaT;

  double _temperature;

  RealGradient _periodicity;
};



inline
Periodic::Periodic(const ModelOptions& options) :
  BoltzmannBoundaryModel(options),
  _deltaT(0.0),
  _periodicity(0)
{
}



inline
Periodic*
Periodic::create(const ModelOptions& options)
{
  return new Periodic(options);
}



inline
PhysicalModelInterface*
Periodic::create_new(void) const
{
  return new Periodic(get_options());
}

#endif // _POISSONDIRICHLET_H_
