// $Id: PoissonNeumann.h 2029 2010-07-21 13:01:58Z maufder $

#ifndef _POISSONNEUMANN_H_
#define _POISSONNEUMANN_H_

#include "PoissonBoundaryModel.h"

#include "point.h"
#include "tensor_value.h"
#include "vector_value.h"



class Elem;


//! The base class for Poisson boundary conditions
class PoissonNeumann : public PoissonBoundaryModel
{

  public:

    //! Destructor
    ~PoissonNeumann(void) {};

    //! Creator function
    static PoissonNeumann* create(const ModelOptions& options);


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
    PoissonNeumann(const ModelOptions& options);

    //! The normal field on the boundary
    double _field;

};



inline
PoissonNeumann::PoissonNeumann(const ModelOptions& options) :
  PoissonBoundaryModel(options),
  _field(0)
{
}



inline
PoissonNeumann*
PoissonNeumann::create(const ModelOptions& options)
{
  return new PoissonNeumann(options);
}



#endif // _POISSONNEUMANN_H_
