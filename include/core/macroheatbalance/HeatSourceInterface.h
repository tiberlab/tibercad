#ifndef _HEATSOURCEINTERFACE_H_
#define _HEATSOURCEINTERFACE_H_

#include "PhysicalModelInterface.h"
#include "SimulationInterface.h"       




//! This class handles the heat source


class HeatSourceInterface : public PhysicalModelInterface
{

public:
  
  virtual ~HeatSourceInterface(void);

  
  //!Get total heat source 
  virtual void get_heat_sources(std::vector<Point> h_point,const Elem* elem,
				std::vector< std::vector<double> >& heat_source){};

  //!Get total heat flux  
  virtual void get_power_fluxes(std::vector<Point> h_point,const Elem* elem,
				std::vector< std::vector<RealGradient>  >& heat_source, bool check_boundary){};
                     

  static HeatSourceInterface* create(const std::string& name,
	       const ModelOptions& options = ModelOptions());

  //! Get the heat source legend
  virtual std::vector<std::string>  get_source_legend(void){};

  //! Get the flux legend
  virtual std::vector<std::string>  get_flux_legend(void){};

	

private:



 
protected:

   HeatSourceInterface(void);
 
 
};


inline
HeatSourceInterface::HeatSourceInterface(void)
{
}

inline
HeatSourceInterface::~HeatSourceInterface(void)
{
}



inline
HeatSourceInterface*
HeatSourceInterface::create(const std::string& name,
				    const ModelOptions& options)
{
  return dynamic_cast<HeatSourceInterface*>(
		PhysicalModelInterface::create(name, options));
}
 


#endif
