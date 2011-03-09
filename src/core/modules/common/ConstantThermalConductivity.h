// $Id: ConstantThermalConductivity.h 2132 2010-10-29 08:19:25Z maufder $

#ifndef _CONSTANTTHERMALCONDUCTIVITY_H_
#define _CONSTANTTHERMALCONDUCTIVITY_H_

#include "ThermalConductivityModel.h"

#include "point.h"
#include "tensor_value.h"
#include "vector_value.h"
#include "PhysicalModelInterface.h"
#include "tiber_dll.h"

class Elem;


//! The base class for Poisson boundary conditions
class TBDLLOCAL ConstantThermalConductivity : public ThermalConductivityModel
{
  
public:
  
  //! Destructor
  virtual ~ConstantThermalConductivity(void) {};
  
  //! Creator function
  static ConstantThermalConductivity* create(const ModelOptions& options);
  
  virtual void calculate(const Elem* elem, const Point& point){};
  
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
  virtual void read_database(void);
  
  
  /* We do not use this here: */
  //  virtual void read_interface_database(void);
  
  
  //! Create a new object of the same type
  virtual PhysicalModelInterface* create_new(void) const;
  
  
private:
  
  RealGradient _kappa;

  //! Constructor
  ConstantThermalConductivity(const ModelOptions& options);
  
};


inline
PhysicalModelInterface*
ConstantThermalConductivity::create_new(void) const
{
  return new   ConstantThermalConductivity(get_options());
}

inline
ConstantThermalConductivity*
ConstantThermalConductivity::create(const ModelOptions& options)
{
  return new  ConstantThermalConductivity(options);
}



#endif // _GRAYMODEL_H_

  
