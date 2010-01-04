// $Id$

#include "HeatModel.h"
#include "Material.h"
#include "LatticeThermalConductivity.h"
#include "SimulationInterface.h"
#include "HeatSourceInterface.h"



HeatModel::HeatModel() :
  kappa(NULL),
  _elem(NULL),
  _lattice_thermal_conductivity(0),
  _heat_source_interface(NULL)
{
  tg = 0.0;
  vg = 0.0;
  cg = 0.0;
}




HeatModel::~HeatModel()
{

  destroy(kappa);
  clear_heat_sources();


}

//==========================================================================//

PhysicalModelInterface* HeatModel::create_new (void) const
{
  return new HeatModel();
}

//==========================================================================//
void HeatModel::create_submodels()
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



  destroy(kappa);

  it = get_options().submodels_begin("Lattice_thermal_conductivity");
  end = get_options().submodels_end("Lattice_thermal_conductivity");


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

   }
    kappa->set_material(get_material());

    kappa->init();
    
    kappa->get_conductivity(_lattice_thermal_conductivity);

      

}


//==========================================================================//
void HeatModel::do_init_alloy(const PhysicalModelInterface *comp_A, const PhysicalModelInterface *comp_B, double xa)
{

  const HeatModel* matA = dynamic_cast< const HeatModel*> (comp_A);
  const HeatModel* matB = dynamic_cast< const HeatModel*> (comp_B);

  destroy(kappa);
  kappa = create_submodel_alloy(matA->kappa, matB->kappa, xa);
  kappa->get_conductivity(_lattice_thermal_conductivity);

  
  _heat_source_models.clear();
  std::vector<ID> source_ids;
  int n = matA->get_heat_source_IDs(source_ids);
  for (int i = 0; i < n; i++)
  {
    ID id = source_ids[i];
    _heat_source_models[id] = create_submodel_alloy(matA->get_heat_source_model(id),
        matB->get_heat_source_model(id), xa);
  }

  // tg =  alloy(matA->tg, matB->tg, xa);
  vg =  alloy(matA->vg, matB->vg, xa);
  cg =  alloy(matA->cg, matB->cg, xa);
  //  kg=  alloy(matA->kg, matB->kg, xa);

  // std::cout<<get_material()->get_name()<<std::endl;

  //  std::cout<< _lattice_thermal_conductivity<<std::endl;
}

void
HeatModel::do_init(void)
{

  Database& db = get_database();

  db.set_section("phonon_velocity/constant");
  vg = db.get("vg", vg);

  //db.set_section("phonon_scattering/constant");
  //tg = db.get("tau_g", tg);
  
  db.set_section("lattice_thermal_capacity/constant");
  cg = db.get("C", cg);

  //kg = _lattice_thermal_conductivity(1,1);
  
  db.set_section("thermal_conductivity/constant");
  kg = db.get("therm_lat_cond_x", kg);

  tg = 3.0 * kg / (vg * vg * cg);  //s

   std::cout<<get_material()->get_name()<<std::endl;

 
   std::cout<<"tg: "<<tg<<std::endl;
   std::cout<<"cg: "<<cg<<std::endl;
   std::cout<<"vg: "<<vg<<std::endl;
   std::cout<<"kg: "<<kg<<std::endl;
   std::cout<<" "<<std::endl;
   //cg = 1;
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
HeatModel::get_total_heat_source(const Elem*  elem, std::vector<Point> h_point,
		  std::vector<double>& total_heat_source)
{

  ID np = h_point.size();

  // for (ID n = 0; n<np; n++)
  //  h_point[n] *= 1e-7;


  total_heat_source.clear();
  total_heat_source.resize(np);

  outer_source_iterator it_outer = _heat_source_models.begin();
  outer_source_iterator end_outer = _heat_source_models.end();
  inner_source_iterator it_inner, end_inner;

  for ( ; it_outer != end_outer; ++it_outer)
  {
    std::vector<std::map<ID,double> > heat_source;
    (it_outer->second)->get_heat_sources(elem,h_point,heat_source);

    for (ID n = 0;  n < np; ++n)
    { 
     
      it_inner =  heat_source[n].begin();
      end_inner = heat_source[n].end();
    
      for ( ; it_inner != end_inner; ++it_inner)
      {
	total_heat_source[n] += it_inner->second;
	
      }
    }
  }

}


void
HeatModel::clear_heat_sources(void)
{

outer_source_iterator it =   _heat_source_models.begin();
outer_source_iterator end =  _heat_source_models.end();

 for ( ; it != end; ++it)
   destroy(it->second);

 _heat_source_models.clear();

}


int
HeatModel::get_heat_source_IDs(std::vector<ID>& ids) const
{
  int n =  _heat_source_models.size();

  ids.resize(n);

  const_outer_source_iterator it =  _heat_source_models.begin();
  const_outer_source_iterator end = _heat_source_models.end();
  int ctr = 0;
  for ( ; it != end; ++it, ctr++)
    ids[ctr] = (it->first);

  return n;


}

void
HeatModel::do_print_info(void)
{
 
  outer_source_iterator it =  _heat_source_models.begin();
  outer_source_iterator end = _heat_source_models.end();
  for ( ; it != end; ++it)
  {
    std::cout<<"Heat sources:"<<std::endl;
    (it->second)->print_info();
  }
   
  
}

void
HeatModel::re_init(void)
{

  // update_lattice_thermal_conductivity();

}

void
HeatModel::update_lattice_thermal_conductivity(void)
{

 kappa->get_conductivity(_lattice_thermal_conductivity);

}



// void
// HeatModel::get_total_power_flux(std::vector<Point> h_point,
// 				std::vector<RealGradient>& total_power_flux)
// {

//   unsigned int np = h_point.size();

//   total_power_flux.clear();
//   total_power_flux.resize(np);

//   outer_source_iterator it_outer = _heat_source_models.begin();
//   outer_source_iterator end_outer = _heat_source_models.end();


//   ID IDtot = 100;
//   std::set<ID> TotalSet;
//   TotalSet.insert(IDtot);


//   for ( ; it_outer != end_outer; it_outer++)
//   {

//     std::vector<std::map<ID,RealGradient> >  partial_power_fluxes;

//     (it_outer->second)->get_power_fluxes(h_point,TotalSet,partial_power_fluxes);

//     for (ID n = 0;  n < np; ++n)
//       total_power_flux[n] += partial_power_fluxes[n][IDtot];



//    }

// }

