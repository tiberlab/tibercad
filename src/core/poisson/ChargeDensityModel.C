#include "ChargeDensityModel.h"
#include "SimulationEnvironment.h"
#include "Material.h"
#include "getpot.h"
#include "Database.h"
//-------------------------------------------------------------------------//

ChargeDensityModel::ChargeDensityModel() :
  _charge_density(0),
  _elem(NULL)
{
}


void ChargeDensityModel::copy_from(const PhysicalModelInterface *rhs)
{
  const ChargeDensityModel* mod = dynamic_cast<const  ChargeDensityModel*> (rhs);

  _charge_density = mod-> _charge_density;
 
}

//-------------------------------------------------------------------------//


void ChargeDensityModel::calculate_VCA (const PhysicalModelInterface *comp_A, 
                                                const PhysicalModelInterface *comp_B, double xa) 
{ 
  const ChargeDensityModel* modA = dynamic_cast<const ChargeDensityModel*>(comp_A);

  const ChargeDensityModel* modB = dynamic_cast<const ChargeDensityModel*>(comp_B);


   alloy(_charge_density,modA->_charge_density, modB->_charge_density,xa);  

  
}

//-------------------------------------------------------------------------//

//--------------------------------------------------------//
void  ChargeDensityModel::read_database(void)
{

  
}

//---------------------------------------------------------//

void ChargeDensityModel::do_init(void)
{


   std::string chd_sim_name = get_options().get_option("charge_density_simulation", "no_dd");

     _chd_sim = SimulationInterface::find_simulation(chd_sim_name);

     // if ((_chd_sim == NULL) || (chd_sim_name.compare("no_dd") == 1))
     //   throw InitFailedException("Unknown " +  chd_sim_name + " simulation");

      
  
}

void  ChargeDensityModel::re_init(void)
{

  if (_chd_sim != NULL)
  {
    SimulationEnvironment& se = _chd_sim->get_environment();
    
    if  (se.contains_element(_elem))
    {

      //  _charge_density =  _chd_sim->ge

    }
    else
    {
      //   _charge_density = 0.0;
    }
  }
  else
  {

    const  Material* mat = get_material();
 
    _charge_density = mat->get_net_doping_density(); 
     
  }
}






//-------------------------------------------------------------------------//
