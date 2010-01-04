#ifndef _CONSTANTHEATSOURCE_H_
#define _CONSTANTHEATSOURCE_H_

#include "PhysicalModelInterface.h"
#include "HeatSourceInterface.h"       




class ConstantHeatSource : public HeatSourceInterface
{

public:
  
  //!Constructor 
  ConstantHeatSource(){}

   //!Destructor
  ~ConstantHeatSource(){}

  //! set the electron conducibility;
  void get_total_heat_source();

  //! Costructor
  static   ConstantHeatSource* create();
 
  
  
  virtual void get_power_fluxes(const Elem*  elem,std::vector<Point> h_point, std::vector<std::map<ID,RealGradient> >& power_fluxes);


  virtual void get_heat_sources(const Elem*  elem,std::vector<Point> h_point, std::vector<std::map<ID, double> >& heat_sources);     


  //! \copydoc HeatSourceInterface::get_source_legend(const std::set<std::string>& variables);
  virtual std::map<ID,std::string> get_source_legend(const std::set<std::string>& variables);

  //!  \copydoc HeatSourceInterface::get_flux_legend(const std::set<std::string>& variables);
  virtual std::map<ID,std::string> get_flux_legend(const std::set<std::string>& variables){};


  //!Set the current element
  virtual void set_heat_model(HeatModel* heat_model){}; 
  
  //! \copydoc PhysicalModel::do_print_info(void)
  virtual void do_print_info(void);

   //! A class that handles the heat source model option
   class HeatSourceParameters
   {
   public:
     //!Constructor
     HeatSourceParameters(void);
    
     //! name of this heat source model
     std::string model_name;
     
     //! name of the drift diffusion siulation
     std::string dd_simul_name;
    
  };

private:

  double _heat_source;

   

  //!Power flux legend
  std::map<ID,std::string> _flux_legend;

  //!Heat Source legend
  std::map<ID,std::string> _source_legend;

  //!Pointer to drift diffusion simulation
  SimulationInterface* _simul;

  //!Pointer to heat model
  // HeatModel* _heat_model;

  //!Heat source variables for drift diffusion
  std::set<ID> ID_set;

  //!Variable map
  std::map<ID,ID> var_map;

 
protected:

    

  virtual void do_init (void);


  virtual void do_init_alloy (const PhysicalModelInterface *comp_A, const PhysicalModelInterface *comp_B, double xa); 

  virtual PhysicalModelInterface* create_new (void) const;
 
 
 
};



inline
PhysicalModelInterface*  ConstantHeatSource::create_new () const
{
  return (new    ConstantHeatSource() ); 
}

inline
ConstantHeatSource*
ConstantHeatSource::create()
{
  return new ConstantHeatSource(); 
}





#endif
