#include "PoissonModel.h"
#include "Material.h"
#include "ChargeDensityModel.h"
 



PoissonModel::PoissonModel(void) :
  _charge_density(0.0),
  _epsilon(0),
  chd_model(NULL),
  dielectric_model(NULL)
{
}
	
 

PoissonModel::~PoissonModel(void)
{

PhysicalModelInterface::destroy(chd_model);

PhysicalModelInterface::destroy(dielectric_model);

}

//==========================================================================//

PhysicalModelInterface* PoissonModel::create_new (void) const
{

   return new PoissonModel;
}

//==========================================================================//
void PoissonModel::do_init()
{
   

   PhysicalModelInterface::destroy(chd_model);


   ModelOptions::const_submodel_iterator it,end;
   it = get_options().submodels_begin("Charge_Density_Model");
   end = get_options().submodels_end("Charge_Density_Model");


   if (it != end)
    {

      chd_model = dynamic_cast<ChargeDensityModel*>(
        PhysicalModelInterface::create("charge_density_model", it->second)); 
  
      if (chd_model == NULL)
	throw InitFailedException("Could not create charge density model");
      
      
    }
   else
   {
    
     chd_model = dynamic_cast<ChargeDensityModel*>(
	 PhysicalModelInterface::create("charge_density_model")); 
     
  
   }

  chd_model->set_material(get_material());

  chd_model->init();

  PhysicalModelInterface::destroy(dielectric_model);

   it = get_options().submodels_begin("Dielectric_Model");
   end = get_options().submodels_end("Dielectric_Model");


   if (it != end)
    {

      dielectric_model = dynamic_cast< DielectricModel*>(
        PhysicalModelInterface::create("dielectric_model", it->second)); 
  
      if ( dielectric_model == NULL)
	throw InitFailedException("Could not create dielectric model");
      
      
    }
   else
   {
     
     dielectric_model = dynamic_cast<DielectricModel*>(
	 PhysicalModelInterface::create("dielectric_model")); 
     
  
   }

   dielectric_model ->set_material(get_material());
   
   dielectric_model->init();
   
   _epsilon = dielectric_model->get_dielectric_constant();
}




//==========================================================================//
void PoissonModel::calculate_VCA (const PhysicalModelInterface *comp_A, const PhysicalModelInterface *comp_B, double xa)
{

  const PoissonModel* matA = dynamic_cast< const PoissonModel*> (comp_A);
  
  const PoissonModel* matB = dynamic_cast< const PoissonModel*> (comp_B);
  
  
  // _epsilon -> build_alloy(matA->get_dielectric, matB->chd_model, xa);

  //   _charge_density -> build_alloy(matA->chd_model, matB->chd_model, xa);
 
}



void  PoissonModel::re_init()
{
  
  update_charge_density();

  //  update_dielectric_constant();      

}

void  PoissonModel::update_charge_density()
{ 

 
  chd_model->set_element(_elem);

  chd_model->re_init();

  _charge_density = chd_model->get_charge_density();

}

//void  PoissonModel::update_dielectric_constant()
//{ 

 
//  chd_model->set_element(_elem);

//  dielectric_model->re_init();

//  charge_density = chd_model->get_charge_density();


//}


