// $Id: LinearThermalConductivity.h 2069 2010-09-08 18:08:39Z gromano $

#ifndef _LINEARTHERMALCONDUCTIVITY_H_
#define _LINEARTHERMALCONDUCTIVITY_H_

#include "point.h"
#include "tensor_value.h"
#include "vector_value.h"

#include "tiber_dll.h"

#include "elem.h"
#include "ThermalConductivityModel.h"






//! The base class for Poisson boundary conditions
class TBDLLOCAL LinearThermalConductivity : public ThermalConductivityModel
{
  
public:
  
  //! Destructor
  virtual ~LinearThermalConductivity(void) {};
  
  //! Creator function
  static LinearThermalConductivity* create(const ModelOptions& options);
  
 virtual void calculate(const libMesh::Elem* elem, const libMesh::Point& point, double temperature);
  
protected:
  
  //! Initialize
  virtual void do_init(void);
  
  /* In some cases it might be useful to reimplement this: */
  // virtual void do_init_interface(const PhysicalModelInterface* comp_A,
  //         const PhysicalModelInterface* comp_B);
  
  virtual void do_init_alloy (const PhysicalModelInterface *comp_A,
			      const PhysicalModelInterface *comp_B, double xa){};
  
  virtual void  read_database_alloy(void){};
  /* This is not used here: */
  virtual void read_database(void){};
  
  
  /* We do not use this here: */
  //  virtual void read_interface_database(void);
  
  
  
private:
   
   double kx0;
   double kz0;
   double mx;
   double mz;
   double z0;
  
  //! Constructor
  LinearThermalConductivity(const ModelOptions& options);
  
};




inline
LinearThermalConductivity*
LinearThermalConductivity::create(const ModelOptions& options)
{
  return new  LinearThermalConductivity(options);
}




#endif // _GRAYMODEL_H_
