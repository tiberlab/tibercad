// $Id: ThermalBoundaryModel.h 2070 2010-09-09 14:51:29Z gromano $

#ifndef _BOLTZMANNBOUNDARYMODEL_H_
#define _BOLTZMANNBOUNDARYMODEL_H_

#include "PhysicalModel.h"


class Elem;
class Point;


//! The base class for thermal balance boundary conditions
class BoltzmannBoundaryModel : public PhysicalModel
{

  public:

    //! Destructor
    ~BoltzmannBoundaryModel(void) {};

    //! Creator function
    static BoltzmannBoundaryModel* create(const ModelOptions& options);


    //! Calculate for a point on the given side
    virtual void calculate(const Elem* elem, unsigned int side,
        const Point& point) = 0;

   void get_coefficients(double& a, double& b, double& c);


  ///!Set the current element
   void set_element(const Elem* elem);
  protected:

    //! Constructor
    BoltzmannBoundaryModel(const ModelOptions& options);

  void set_coefficients(double a, double b, double c);


  private:

  double _alpha;
  double _beta;
  double _gamma;

};


inline
BoltzmannBoundaryModel::BoltzmannBoundaryModel(const ModelOptions& options) :
  PhysicalModel(options)
{
}


inline
void
BoltzmannBoundaryModel::get_coefficients(double& a, double& b, double& c)
{
 a = _alpha;
 b = _beta;
 c = _gamma;
}


inline
void
BoltzmannBoundaryModel::set_coefficients(double a, double b, double c)
{
 _alpha = a;
 _beta  = b;
 _gamma = c;
}


#endif // _THERMALBOUNDARYMODEL_H_
