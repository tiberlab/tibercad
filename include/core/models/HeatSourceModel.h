// $Id: HeatSourceModel.h 2457 2011-03-06 23:52:12Z gromano $

#ifndef _HEATSOURCEMODEL_H_
#define _HEATSOURCEMODEL_H_

#include "PhysicalModelInterface.h"

#include "tensor_value.h"

//#undef  TIBER_MODULE_PREFIX
//#define TIBER_MODULE_PREFIX thermal_conductivity

class Elem;
class Point;

using namespace std;

//! The base class for Poisson boundary conditions
class HeatSourceModel : public PhysicalModelInterface
{

  public:

    //! Destructor
    ~HeatSourceModel(void) {};

     //! Creator function
   static HeatSourceModel* create(const ModelOptions& options);

  Real get_heat_source(void) const;

  virtual void calculate(const Elem* elem, const Point& point){};

  protected:

    //! Constructor
  HeatSourceModel(const ModelOptions& options);

//! Calculate for a point on the given side
  //void calculate(const Elem* elem, const Point& point);

  void set_heat_source(Real heat_source);

  private:

  Real _heat_source;

};

inline
Real
HeatSourceModel::get_heat_source(void) const
{

  return _heat_source;
}


inline 
void 
HeatSourceModel::set_heat_source(Real heat_source)
{
  _heat_source = heat_source;
}



inline
HeatSourceModel::HeatSourceModel(const ModelOptions& options) :
PhysicalModelInterface(options)
{
}




#endif // _THERMALCONDUCTIVITYMODEL_H_
