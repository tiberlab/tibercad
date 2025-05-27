// $Id$

#ifndef _WIBOUNDARYMODEL_H_
#define _WIBOUNDARYMODEL_H_

#include "PhysicalModel.h"


class Elem;
class Point;


//! The base class for WI boundary conditions
class WIBoundaryModel : public PhysicalModel
{

  public:

    //! Destructor
    ~WIBoundaryModel(void) {};

    //! Creator function
    static WIBoundaryModel* create(const MaterialBoundary* boundary,
        const ModelOptions& options);


    //! Calculate for a point on the given side
    virtual void calculate(const Elem* elem, unsigned int side,
        const Point& point) = 0;

    void get_coefficients(double& a, double& b, double& c);


  protected:

    //! Constructor
    WIBoundaryModel(const ModelOptions& options);

    void set_coefficients(double a, double b, double c);


  private:

    double _alpha;
    double _beta;
    double _gamma;

};



inline
WIBoundaryModel::WIBoundaryModel(const ModelOptions& options) :
  PhysicalModel(options),
  _alpha(0),
  _beta(1),
  _gamma(0)
{
}



inline
void
WIBoundaryModel::get_coefficients(double& a, double& b, double& c)
{
  a = _alpha;
  b = _beta;
  c = _gamma;
}


inline
void
WIBoundaryModel::set_coefficients(double a, double b, double c)
{
  _alpha = a;
  _beta = b;
  _gamma = c;
}


#endif // _WIBOUNDARYMODEL_H_
