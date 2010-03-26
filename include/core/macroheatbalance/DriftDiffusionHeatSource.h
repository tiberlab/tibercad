#ifndef _DRIFTDIFFUSIONHEATSOURCE_H_
#define _DRIFTDIFFUSIONHEATSOURCE_H_

#include "PhysicalModelInterface.h"
#include "HeatSourceInterface.h"




class DriftDiffusionHeatSource : public HeatSourceInterface
{

public:

   //!Destructor
  ~DriftDiffusionHeatSource(){};

  //! set the electron conducibility;
  void get_total_heat_source();

  //! Costructor
  static  DriftDiffusionHeatSource* create(const ModelOptions& options);



  virtual void get_power_fluxes(const Elem*  elem, std::vector<Point> h_point, std::vector<std::map<ID,RealGradient> >& power_fluxes);


  virtual void get_heat_sources( const Elem*  elem, std::vector<Point> h_point, std::vector<std::map<ID, double> >& heat_sources);


  //! \copydoc HeatSourceInterface::get_source_legend(const std::set<std::string>& variables);
  virtual std::map<ID,std::string> get_source_legend(const std::set<std::string>& variables);

  //!  \copydoc HeatSourceInterface::get_flux_legend(const std::set<std::string>& variables);
  virtual std::map<ID,std::string> get_flux_legend(const std::set<std::string>& variables);


  //!Set the current element
  virtual void set_heat_model(HeatModel* heat_model);

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

  //!Heat model option
  HeatSourceParameters heat_source_opt;


  enum heat_variables
    {
      EJOULE = 0,
      HJOULE,
      RECHEAT,
      EPELTH,
      HPELTH,
      WNX,
      WPX
    };

  //!Power flux legend
  std::map<ID,std::string> _flux_legend;

  //!Heat Source legend
  std::map<ID,std::string> _source_legend;

  //!Pointer to drift diffusion simulation
  SimulationInterface* _simul;

  //!Pointer to heat model
  HeatModel* _heat_model;

  //!Heat source variables for drift diffusion
  std::set<ID> ID_set;

  //!Variable map
  std::map<ID,ID> var_map;


protected:


  //!Constructor
  DriftDiffusionHeatSource(const ModelOptions& options);


  virtual void do_init (void);


  virtual PhysicalModelInterface* create_new (void) const;



};



inline
DriftDiffusionHeatSource::DriftDiffusionHeatSource(const ModelOptions& options) :
  HeatSourceInterface(options)
{
}



inline
PhysicalModelInterface*  DriftDiffusionHeatSource::create_new () const
{
  return (new    DriftDiffusionHeatSource(get_options()) );
}

inline
DriftDiffusionHeatSource*
DriftDiffusionHeatSource::create(const ModelOptions& options)
{
  return new DriftDiffusionHeatSource(options);
}





inline
void
DriftDiffusionHeatSource::set_heat_model(HeatModel* heat_model)
{

 _heat_model = heat_model;

}

// inline
// std::map<ID,std::string>
// DriftDiffusionHeatSource::get_source_legend(void)
// {

//  return  _source_legend;

// }

// inline
// std::map<ID,std::string>
// DriftDiffusionHeatSource::get_flux_legend(void)
// {


//  return  _flux_legend;

// }


#endif
