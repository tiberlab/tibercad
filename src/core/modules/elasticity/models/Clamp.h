// $Id$

#ifndef _CLAMP_H_
#define _CLAMP_H_

#include "BoundaryModel.h"

#include "point.h"
#include "tensor_value.h"
#include "vector_value.h"



class Elem;


//! The base class for Poisson boundary conditions
class Clamp : public BoundaryModel
{

  public:

    //! Destructor
    ~Clamp(void) {};

    //! Creator function
    static Clamp* create(const ModelOptions& options);


    //! Calculate for a point on the given side
    virtual void calculate(const Elem* elem, unsigned int side,
			   const Point& point){};


  protected:

    //! Initialize
  virtual void do_init(void);


  private:

    //! Constructor
    Clamp(const ModelOptions& options);

    //! The boundary potential
    //double _potential;

};



inline
Clamp::Clamp(const ModelOptions& options) :
  BoundaryModel(options)
{
}



inline
Clamp*
Clamp::create(const ModelOptions& options)
{
  return new Clamp(options);
}




#endif // _POISSONDIRICHLET_H_
