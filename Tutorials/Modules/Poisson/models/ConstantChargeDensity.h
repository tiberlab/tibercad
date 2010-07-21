// $Id$

#ifndef _CONSTANTCHARGEDENSITY_H_
#define _CONSTANTCHARGEDENSITY_H_

#include "MyChargeDensityModel.h"


class Elem;
class Point;

// Base class for charge density models
class ConstantChargeDensity : public ChargeDensityModel
{

  public:

    virtual ~ConstantChargeDensity(void) {};

    static ConstantChargeDensity* create(const ModelOptions& options);


  protected:

    ConstantChargeDensity(const ModelOptions& options);

    virtual void do_init(void);

    virtual double calculate_charge_density(const Elem* elem, const Point& point);

  private:

    double _charge_density;

};


inline
ConstantChargeDensity::ConstantChargeDensity(const ModelOptions& options) :
  ChargeDensityModel(options),
  _charge_density(0.0)
{
}


inline
ConstantChargeDensity*
ConstantChargeDensity::create(const ModelOptions& options)
{
  return new ConstantChargeDensity(options);
}


#endif // _CHARGEDENSITYMODEL_H_
