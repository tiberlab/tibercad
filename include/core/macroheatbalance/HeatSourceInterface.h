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
  virtual void get_heat_source(std::vector<Point> h_point,const Elem* elem,
			       std::vector< double >& heat_source){};

  //!Get total heat flux  
  virtual void get_flux_heat_source(std::vector<Point> h_point,const Elem* elem,
			       std::vector< RealGradient  >& heat_source){};

  //!Get heat source contributions  
  virtual void get_heat_source_output(std::vector<Point> h_point,const Elem* elem,
				     std::vector< std::vector< double > >& heat_source){};

  
  static HeatSourceInterface* create(const std::string& name,
	       const ModelOptions& options = ModelOptions());

  //! Get the legend
  virtual std::vector<std::string>  get_legend(void){};
	

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
