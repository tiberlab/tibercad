// $Id: ConstantChargeDensity.h 3864 2014-06-13 10:24:26Z maufder $

#ifndef _CONSTANTCHARGEDENSITY_H_
#define _CONSTANTCHARGEDENSITY_H_

#include "ChargeDensityModel.h"
#include "SolutionProvider.h"

class Elem;
class Point;

// Base class for charge density models
class ChargeDensity : public ChargeDensityModel
{

  public:

    virtual ~ChargeDensity(void) {};

    static ChargeDensity* create(const ModelOptions& options);


  protected:

    ChargeDensity(const ModelOptions& options);

    virtual void do_init(void);

    virtual double calculate_charge_density(const Elem* elem, const Point& point);

  private:

    double _charge_density;
   
    std::vector<SolutionProvider> _charge_density_sim;
};


inline
ChargeDensity::ChargeDensity(const ModelOptions& options) :
  ChargeDensityModel(options),
  _charge_density(0.0)
{
}


inline
ChargeDensity*
ChargeDensity::create(const ModelOptions& options)
{
  ChargeDensity* cd = new ChargeDensity(options);
  std::cout << "Charge Density Created ID: "<< cd->get_id()<< std::endl;

  return cd;

}



#endif // _CHARGEDENSITYMODEL_H_
