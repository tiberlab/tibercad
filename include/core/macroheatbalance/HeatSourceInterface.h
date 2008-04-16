#ifndef _HEATSOURCEINTERFACE_H_
#define _HEATSOURCEINTERFACE_H_

#include "PhysicalModelInterface.h"
#include "SimulationInterface.h"       
//#include "HeatModel.h" 

class HeatModel;

//! This class handles the heat source


class HeatSourceInterface : public PhysicalModelInterface
{

public:
  
  virtual ~HeatSourceInterface(void);

  
  //!Get total heat source 
  virtual void get_heat_sources(std::vector<Point> h_point, const std::set<ID>& ids,
				std::vector<std::map<ID, double> >& heat_sources){};

  //!Get total heat flux  
  virtual void get_power_fluxes(std::vector<Point> h_point, const std::set<ID>& ids,
				std::vector<std::map<ID,RealGradient> >& heat_source){};
                     

  static HeatSourceInterface* create(const std::string& name,
	       const ModelOptions& options = ModelOptions());

  virtual std::map<ID,std::string> get_source_legend(const std::set<std::string>& variables){};

  //! Get the flux legend
  virtual std::map<ID,std::string>  get_flux_legend(const std::set<std::string>& variables){};


  //! Set the heat model
  virtual void set_heat_model(HeatModel* _heat_model){};



  

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
