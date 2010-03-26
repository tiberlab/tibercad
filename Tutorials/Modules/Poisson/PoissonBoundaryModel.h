// $Id: PoissonBoundaryModel.h 1856 2010-03-22 15:55:26Z maufder $

#ifndef _POISSONBOUNDARYMODEL_H_
#define _POISSONBOUNDARYMODEL_H_

#include "PhysicalModel.h"



#undef TIBER_MODULE_PREFIX
#define TIBER_MODULE_PREFIX mypoisson

class Elem;
class Point;


//! The base class for Poisson boundary conditions
class PoissonBoundaryModel : public PhysicalModel
{

  public:

    //! Destructor
    ~PoissonBoundaryModel(void) {};

    //! Creator function
    static PoissonBoundaryModel* create(const ModelOptions& options);


    //! Calculate for a point on the given side
    virtual void calculate(const Elem* elem, unsigned int side,
        const Point& point) = 0;

    void get_coefficients(double& a, double& b, double& c);


  protected:

    //! Constructor
    PoissonBoundaryModel(const ModelOptions& options);

    void set_coefficients(double a, double b, double c);


  private:

    double _alpha;
    double _beta;
    double _gamma;

};



inline
PoissonBoundaryModel::PoissonBoundaryModel(const ModelOptions& options) :
  PhysicalModel(options),
  _alpha(0),
  _beta(1),
  _gamma(0)
{
}



inline
void
PoissonBoundaryModel::get_coefficients(double& a, double& b, double& c)
{
  a = _alpha;
  b = _beta;
  c = _gamma;
}


inline
void
PoissonBoundaryModel::set_coefficients(double a, double b, double c)
{
  _alpha = a;
  _beta = b;
  _gamma = c;
}


#endif // _POISSONBOUNDARYMODEL_H_
