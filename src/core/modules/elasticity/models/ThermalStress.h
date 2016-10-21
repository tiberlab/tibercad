// $Id: ThermalStress.h 2451 2011-03-05 14:45:46Z maufder $

#ifndef _THERMALSTRESS_H_
#define _THERMALSTRESS_H_

#include "BodyForceModel.h"

#include "TemperatureInterface.h"
#include "vector_value.h"
#include "tiber_dll.h"

class Elem;



//! The base class for Poisson boundary conditions
class TBDLLOCAL ThermalStress : public BodyForceModel
{

  public:
 
  //! Destructor
  ~ThermalStress(void);
  
  //! Creator function
  static ThermalStress* create(const ModelOptions& options);
  
  virtual void calculate(const libMesh::Elem* elem, const libMesh::Point& point);

  protected:

    //! Initialize
    virtual void do_init(void);


    //! Read expansion coefficients
    virtual void read_database(void);


  private:
  

    //! Constructor
    ThermalStress(const ModelOptions& options);

    //! Thermal expansion coefficients for the crystal directions
    libMesh::RealVectorValue _alpha;

    //! A reference temperature
    double _ref_temp;

    //! From where to get the temperature
    TemperatureInterface _temp;


  
};




inline
ThermalStress*
ThermalStress::create(const ModelOptions& options)
{ 
  return new  ThermalStress(options);
}




#endif // _THERMALSTRESS_H_
