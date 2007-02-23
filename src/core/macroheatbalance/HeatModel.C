#include "HeatModel.h"
#include "Material.h"
HeatModel::HeatModel() 
{
  kappa = NULL;  

};


HeatModel::~HeatModel()
{

  PhysicalModelInterface::destroy(kappa);

}

//==========================================================================//

PhysicalModelInterface* HeatModel::create_new (void) const
{
  return new HeatModel();
}

//==========================================================================//
void HeatModel::do_init()
{

  PhysicalModelInterface::destroy(kappa);
  
  const ModelOptions& opt =  get_options ();

  kappa =dynamic_cast<LatticeThermalConductivity*>
    (  PhysicalModelInterface::create("lat_therm_cond_" + get_material()->get_structure(),  opt  )  );
  
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

