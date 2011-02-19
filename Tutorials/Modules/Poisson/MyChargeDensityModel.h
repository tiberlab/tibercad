// $Id$

#ifndef _CHARGEDENSITYMODEL_H_
#define _CHARGEDENSITYMODEL_H_

#include "PhysicalModelInterface.h"


class Elem;
class Point;

// Base class for charge density models
class MyChargeDensityModel : public PhysicalModelInterface
{

  public:

    virtual ~MyChargeDensityModel(void) {};

    double get_charge_density(const Elem* elem, const Point&);


  protected:

    MyChargeDensityModel(const ModelOptions& options);

    virtual double calculate_charge_density(const Elem* elem, const Point& point) = 0;

};


MyChargeDensityModel::MyChargeDensityModel(const ModelOptions& options) :
  PhysicalModelInterface(options)
{
}


double
MyChargeDensityModel::get_charge_density(const Elem* elem, const Point& point)
{
  return calculate_charge_density(elem, point);
}


#endif // _CHARGEDENSITYMODEL_H_
