// $Id: ChargeDensityModel.h 3867 2014-06-16 06:54:04Z maufder $

#ifndef _CHARGEDENSITYMODEL_H_
#define _CHARGEDENSITYMODEL_H_

#include "PhysicalModelInterface.h"


class Elem;
class Point;

// Base class for charge density models
class ChargeDensityModel : public PhysicalModelInterface
{

  public:

    virtual ~ChargeDensityModel(void) {};

    double get_charge_density(const Elem* elem, const Point&);


  protected:

    ChargeDensityModel(const ModelOptions& options);

    virtual double calculate_charge_density(const Elem* elem, const Point& point) = 0;

};


ChargeDensityModel::ChargeDensityModel(const ModelOptions& options) :
  PhysicalModelInterface(options)
{
}


double
ChargeDensityModel::get_charge_density(const Elem* elem, const Point& point)
{
  return calculate_charge_density(elem, point);
}


#endif // _CHARGEDENSITYMODEL_H_
