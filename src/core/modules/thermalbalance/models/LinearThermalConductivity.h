// $Id$

#ifndef _LINEARTHERMALCONDUCTIVITY_H_
#define _LINEARTHERMALCONDUCTIVITY_H_

#include "ThermalConductivityModel.h"

#include "point.h"
#include "tensor_value.h"
#include "vector_value.h"
#include "PhysicalModelInterface.h"


class Elem;


//! The base class for Poisson boundary conditions
class LinearThermalConductivity : public ThermalConductivityModel
{
  
public:
  
  //! Destructor
  virtual ~LinearThermalConductivity(void) {};
  
  //! Creator function
  static LinearThermalConductivity* create(const ModelOptions& options);
  
  virtual void calculate(const Elem* elem, const Point& point);
  
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
  
  
  //! Create a new object of the same type
  virtual PhysicalModelInterface* create_new(void) const;
  
  
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
PhysicalModelInterface*
LinearThermalConductivity::create_new(void) const
{
  return new   LinearThermalConductivity(get_options());
}

inline
LinearThermalConductivity*
LinearThermalConductivity::create(const ModelOptions& options)
{
  return new  LinearThermalConductivity(options);
}




#endif // _GRAYMODEL_H_
