#include "HeatModel.h"
#include "Material.h"
#include "LatticeThermalConductivity.h"
#include "SimulationInterface.h"
#include "HeatSourceInterface.h"



HeatModel::HeatModel() :
  kappa(NULL),
  kappa_carrier(NULL),
  _elem(NULL),
  _lattice_thermal_conductivity(0),
  _electrons_thermal_conductivity(0),
  _heat_source_interface(NULL),
  _holes_thermal_conductivity(0)
{
}
	
 
 

HeatModel::~HeatModel()
{
  
  PhysicalModelInterface::destroy(kappa);
  PhysicalModelInterface::destroy(kappa_carrier);
  


  clear_heat_sources();
  //clear_thermal_conductivity();

}

//==========================================================================//

PhysicalModelInterface* HeatModel::create_new (void) const
{
  return new HeatModel();
}

//==========================================================================//
void HeatModel::do_init()
{

  ModelOptions::const_submodel_iterator it;
  ModelOptions::const_submodel_iterator end;
  
  //Heat source models
  it = get_options().submodels_begin("heat_source");
  end = get_options().submodels_end("heat_source");

  for ( ; it != end; ++it)
  {
    const std::string& name = (it->second).get_option("model", "");
    add_heat_source_model(name,it->second);
  }

//   //Thermal conductivity models
//   it = get_options().submodels_begin("thermal_conductivity");
//   end = get_options().submodels_end("thermal_conductivity");
  
//   for ( ; it != end; ++it)
//   {
//     const std::string& name = (it->second).get_option("type", "");
//     add_thermal_conductivity_model(name,it->second);
//   }





  PhysicalModelInterface::destroy(kappa);

  it = get_options().submodels_begin("Lattice_thermal_condictivity");
  end = get_options().submodels_end("Lattice_thermal_condictivity");


  if (it != end)
  {
    
    
    kappa =dynamic_cast<LatticeThermalConductivity*>(
		      PhysicalModelInterface::create("lat_therm_cond_" +
		      get_material()->get_structure(), it->second)); 
    
  
   if (kappa == NULL)
      throw InitFailedException("Could not create lattice thermal conductivity model");

 
    }
   else
   {
     kappa = dynamic_cast<LatticeThermalConductivity*>(
	  PhysicalModelInterface::create("lat_therm_cond_" +
		get_material()->get_structure()));
  
     //std::cout<<get_material()->get_structure()<<std::endl;

    }
 
  // std::cout<<get_material()->get_structure()<<std::endl;

    kappa->set_temperature(SimulationOptions::temperature);
    
    kappa->set_material(get_material());

    kappa->init();




   
}



//==========================================================================//
void HeatModel::calculate_VCA (const PhysicalModelInterface *comp_A, const PhysicalModelInterface *comp_B, double xa)
{

  const HeatModel* matA = dynamic_cast< const HeatModel*> (comp_A);

  const HeatModel* matB = dynamic_cast< const HeatModel*> (comp_B);

  kappa->build_alloy(matA->kappa, matB->kappa, xa);

  

}



void HeatModel::re_init()
{

  update_lattice_thermal_conductivity();
      
}



void HeatModel::update_lattice_thermal_conductivity()
{
  //Insert phase
  kappa->set_temperature(_temperature);

  //Update phase
  kappa->update_tensor();

  //Getting result phase
  kappa->get_conductivity(_lattice_thermal_conductivity);

}


void
HeatModel::add_heat_source_model(const std::string& model_name, 
                                const ModelOptions& options)
{


  HeatSourceInterface* model =
      HeatSourceInterface::create(model_name, options);
 
  if (model == NULL)
    throw InitFailedException("No such heat source model: " + model_name);
 
  ID id = model->get_id();
  _heat_source_models[id] = model;

  model->set_heat_model(this);
  model->set_material(get_material());
  model->set_simulator_id(get_simulator_id());
  model->init();


}

void
HeatModel::add_thermal_conductivity_model(const std::string& model_name, 
                                const ModelOptions& options)
{

  //  std::cout<< get_material()->get_structure()<<std::endl;
  // if (model_name == "lattice")
  //{   model_name = "lat_therm_cond_" + get_material()->get_structure();}
	  

  ThermalConductivityInterface* model =
      ThermalConductivityInterface::create(model_name, options);
 
  if (model == NULL)
    throw InitFailedException("No such thermal conductivity model: " + model_name);
 
  ID id = model->get_id();
  _thermal_conductivity_models[id] = model;
  model->set_material(get_material());
  model->set_simulator_id(get_simulator_id());
  model->init();



}



void
HeatModel::get_total_heat_source(std::vector<Point> h_point,
		  std::vector<double>& total_heat_source)
{

  ID np = h_point.size();

  total_heat_source.clear();
  total_heat_source.resize(np);
 
  outer_source_iterator it_outer = _heat_source_models.begin();
  outer_source_iterator end_outer = _heat_source_models.end();
  
  ID IDtot = 100;
  std::set<ID> TotalSet;
  TotalSet.insert(IDtot);

  for ( ; it_outer != end_outer; ++it_outer)
  {

    //std::vector<std::vector<double> >  partial_heat_source;

    std::vector<std::map<ID,double> >  partial_heat_source;

    (it_outer->second)->get_heat_sources(h_point,TotalSet,partial_heat_source);

    for (ID n = 0;  n < np; ++n)
      total_heat_source[n] = partial_heat_source[n][IDtot];


  }
}



void 
HeatModel::get_total_power_flux(std::vector<Point> h_point,
				std::vector<RealGradient>& total_power_flux)
{
  
  unsigned int np = h_point.size();

  //
  total_power_flux.resize(np);
 
  outer_source_iterator it_outer = _heat_source_models.begin();
  outer_source_iterator end_outer = _heat_source_models.end();


  ID IDtot = 100;
  std::set<ID> TotalSet;
  TotalSet.insert(IDtot);
  
  for ( ; it_outer != end_outer; it_outer++)
  {

    std::vector<std::map<ID,RealGradient> >  partial_power_fluxes;
    
    (it_outer->second)->get_power_fluxes(h_point,TotalSet,partial_power_fluxes);

    for (ID n = 0;  n < np; ++n) 
      total_power_flux[n] = partial_power_fluxes[n][IDtot];
    
      
   }
 
}


void
HeatModel::clear_heat_sources(void)
{

outer_source_iterator it =   _heat_source_models.begin();
outer_source_iterator end =  _heat_source_models.end();

for ( ; it != end; ++it)
 {
  
     PhysicalModelInterface::destroy(it->second);

 }

     _heat_source_models.clear();

}

void
HeatModel::clear_thermal_conductivity(void)
{
                                 

outer_conductivity_iterator it =  _thermal_conductivity_models.begin();
outer_conductivity_iterator end = _thermal_conductivity_models.end();

for ( ; it != end; ++it)
 {
  
     PhysicalModelInterface::destroy(it->second);

 }

     _thermal_conductivity_models.clear();

}


int
HeatModel::get_heat_source_IDs(std::vector<ID>& ids)
{
  int n =  _heat_source_models.size();

  ids.resize(n);

  outer_source_iterator it =  _heat_source_models.begin();
  outer_source_iterator end = _heat_source_models.end();
  int ctr = 0;
  for ( ; it != end; ++it, ctr++)
    ids[ctr] = (it->first);

  return n;


}

void
HeatModel::do_print_info(void)
{
  std::string space = "        ";

  std::cout<<space<<"Heat sources:"<<std::endl;
  outer_source_iterator it =  _heat_source_models.begin();
  outer_source_iterator end = _heat_source_models.end();
  for ( ; it != end; ++it)
    (it->second)->print_info();


}


// void
// HeatModel:: update_heat_source_model()
// {

//   outer_source_iterator it =  _heat_source_models.begin();
//   outer_source_iterator end = _heat_source_models.end();
  
//   for ( ; it != end; ++it, it++ )
//   {

//     //    it->set_element();
//     //it->set_side();

//     //it->re_init();

//   }

// }
