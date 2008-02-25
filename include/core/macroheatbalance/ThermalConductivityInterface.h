#ifndef _THERMALCONDUCTIVITYINTERFACE_H_
#define _THERMALCONDUCTIVITYINTERFACE_H_

#include "PhysicalModelInterface.h"
#include "SimulationInterface.h"       



//! This class handles the heat source


class ThermalConductivityInterface : public PhysicalModelInterface
{

public:
  
  virtual ~ThermalConductivityInterface(void);

  

  virtual void get_thermal_conductivity(std::vector<Point> h_point,const Elem* elem,
			       std::vector< double >& thermal_conductivity){};

  

  virtual void get_thermal_conductivity(std::vector<Point> h_point,const Elem* elem,
				     std::vector< std::vector< double > >& thermal_conductivity){};

  
  static ThermalConductivityInterface* create(const std::string& name,
	       const ModelOptions& options = ModelOptions());

  virtual std::vector<std::string>  get_legend(void){};
	

private:



 
protected:

   ThermalConductivityInterface(void);
 
 
};


inline
ThermalConductivityInterface::ThermalConductivityInterface(void)
{
}

inline
ThermalConductivityInterface::~ThermalConductivityInterface(void)
{
}



inline
ThermalConductivityInterface*
ThermalConductivityInterface::create(const std::string& name,
				    const ModelOptions& options)
{


  return dynamic_cast<ThermalConductivityInterface*>(
		PhysicalModelInterface::create(name,options));
}
 
 

#endif
