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
  virtual void get_heat_sources(std::vector<Point> h_point,
				std::vector< std::vector<double> >& heat_source){};

  //!Get total heat flux  
  virtual void get_power_fluxes(std::vector<Point> h_point,
				std::vector< std::vector<RealGradient>  >& heat_source){};
                     

  static HeatSourceInterface* create(const std::string& name,
	       const ModelOptions& options = ModelOptions());

  //! Get the heat source legend
  virtual std::vector<std::string>  get_source_legend(const std::set<std::string>& variables){};

  //! Get the flux legend
  virtual std::vector<std::string>  get_flux_legend(void){};

   //!Init the heat source model
  //virtual void re_init(void){};

  //!Set the current element
  // virtual void set_coordinate(local_info* coord){};

  virtual void set_heat_model(HeatModel* _heat_model){};

  ////!Set the current element
  // virtual void set_element(const Elem* elem){};
  
  ////!Set the current elemement side index
  //  virtual void set_side(unsigned int side){};	


  

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
