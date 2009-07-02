#include "PoissonModel.h"
#include "Material.h"
#include "ChargeDensityModel.h"
#include "Database.h"
#include "getpot.h"
#include "PyroPolarization.h"


PoissonModel::PoissonModel(void) :
  _pyropolarization(NULL),
  _strain_sim(NULL),
  _polarization(3, 0.0),
  _chd_sim(NULL),
  _epsilon_model(NULL),
  _epsilon(0)
{
}



PoissonModel::~PoissonModel(void)
{

  //destroy(chd_model);
  destroy(_pyropolarization);
  destroy(_epsilon_model);

}

//==========================================================================//

PhysicalModelInterface* PoissonModel::create_new (void) const
{
   return new PoissonModel;
}

//==========================================================================//
void PoissonModel::do_init()
{

  //Pyropolarization



   //Density Charge Model

      model_opt.chd_sim_name = get_options().get_option("charge_density_simulation", "no_sim");

      if  (model_opt.chd_sim_name.compare("no_sim")==1)
      {
	_chd_sim =  SimulationInterface::find_simulation( model_opt.chd_sim_name);

	if (_chd_sim == NULL)
	  throw InitFailedException("Unknown charge_density simulation: " +   model_opt.chd_sim_name );
	else
	  charge_id = _chd_sim->get_variable_id("charge_density");
      }
      //doping

      model_opt.add_doping = get_options().get_option("add_doping",false);
      //---------------





      //  ModelOptions::const_submodel_iterator it,end;

    //    destroy(chd_model); */



/*    it = get_options().submodels_begin("Charge_Density_Model"); */
/*    end = get_options().submodels_end("Charge_Density_Model"); */


/*    if (it != end) */
/*     { */

/*       chd_model = dynamic_cast<ChargeDensityModel*>( */
/*         PhysicalModelInterface::create("charge_density_model", it->second));  */

/*       if (chd_model == NULL) */
/* 	throw InitFailedException("Could not create charge density model"); */


/*     } */
/*    else */
/*    { */

/*      //The Simulation name is given from Physics Section (temporany solution) */

/*      chd_model = dynamic_cast<ChargeDensityModel*>( */
/*      PhysicalModelInterface::create("charge_density_model"),get_options());  */


/*    } */

/*   chd_model->set_material(get_material()); */

/*   chd_model->init(); */


   //Dielectric Model
   //   destroy(dielectric_model);

   // it = get_options().submodels_begin("Dielectric_Model");
   //  end = get_options().submodels_end("Dielectric_Model");


 //   if (it != end)
//     {


//       dielectric_model = dynamic_cast< DielectricModel*>(
//         PhysicalModelInterface::create("dielectric_model", it->second));

//       if ( dielectric_model == NULL)
// 	throw InitFailedException("Could not create dielectric model");


//     }
//    else
//    {

//      dielectric_model = dynamic_cast<DielectricModel*>(
// 	 PhysicalModelInterface::create("dielectric_model"));


//    }

//    dielectric_model ->set_material(get_material());

//    dielectric_model->init();

//   _epsilon = dielectric_model->get_dielectric_constant();



}


void
PoissonModel::create_submodels(void)
{
  model_opt.piezo_pol = get_options().get_option("piezoelectric_field", false);
  if (model_opt.pyro_pol)
  {
    std::string strain_simul = get_options().get_option("strain_simulation", "");
    _strain_if.set_simulation(strain_simul);
  }

  model_opt.pyro_pol = get_options().get_option("piroelecric_field", false);
  if (model_opt.pyro_pol)
  {
    destroy(_pyropolarization);
    _pyropolarization = PyroPolarization::create(get_material());
    _pyropolarization->set_material(get_material());
    _pyropolarization->set_simulator_id(get_simulator_id());
    _pyropolarization->init();
  }


  //Dielectric constant
  destroy(_epsilon_model);

  const ModelOptions& opt =  get_options ();

  _epsilon_model = OptDielectricConstant::create(get_material()->get_structure(), opt);

  _epsilon_model->set_material(get_material());

  _epsilon_model->init();

  _epsilon_model->get_dielectric_real(_epsilon);
}



//==========================================================================//
void PoissonModel::do_init_alloy (const PhysicalModelInterface *comp_A, const PhysicalModelInterface *comp_B, double xa)
{

  //const PoissonModel* matA = dynamic_cast< const PoissonModel*> (comp_A);
  // const PoissonModel* matB = dynamic_cast< const PoissonModel*> (comp_B);

  //destroy(dielectric_model);
  // dielectric_model = create_submodel_alloy(matA->dielectric_model,
  //    matB->dielectric_model,xa);

   //   alloy( _pyropolarization(1),matA->_pyropolarization(1),matB->_pyropolarization(1),xa);
   //alloy( _pyropolarization(2),matA->_pyropolarization(2),matB->_pyropolarization(2),xa);
   //alloy( _pyropolarization(3),matA->_pyropolarization(3),matB->_pyropolarization(3),xa);

}



void  PoissonModel::re_init()
{

  // update_charge_density();

  // update_built_in_polarization();

  //  update_dielectric_constant();

}

void  PoissonModel::update_charge_density()
{


  // chd_model->set_element(_elem);

  // chd_model->re_init();

  //_charge_density = chd_model->get_charge_density();

}

//void  PoissonModel::update_dielectric_constant()
//{


//  chd_model->set_element(_elem);

//  dielectric_model->re_init();

//  charge_density = chd_model->get_charge_density();


//}

//void  PoissonModel::update_built_in_polarization()
//{


// if (model_opt.piezo_pol)
// {

// }

//}


void
PoissonModel::get_charge_density(const std::vector<Point> q_point, std::vector<double>& charge_density)
{


  charge_density.resize(q_point.size() );
  charge_density.clear();

  if (_chd_sim!=NULL)
    if (_chd_sim->get_solution(_elem,q_point,charge_id, charge_density)){}


  if (model_opt.add_doping)
  {
    const Material* mat = get_material();
    for(unsigned int n=0;n<q_point.size();n++)
      charge_density[n] += mat->get_net_doping_density()*Constants::e;
  }


}


const RealVectorValue&
PoissonModel::get_total_polarization(void)
{

  _polarization(0)=0.0;
  _polarization(1)=0.0;
  _polarization(2)=0.0;

  if (model_opt.piezo_pol)
  {
    _strain_if.get_strain_data(_elem, _strain, _polarization);
  }

  if (model_opt.pyro_pol)
  {
     _coord = _elem->centroid();

     double temp = SimulationOptions::temperature;

    _pyropolarization->calculate_polarization(_elem,_coord,temp);
    _polarization(0) += _pyropolarization->get_polarization()(1);
    _polarization(1) += _pyropolarization->get_polarization()(2);
    _polarization(2) += _pyropolarization->get_polarization()(3);
  }

  //Scaling from Q/m2 tp Q/cm2
  _polarization(0) *=1e-4;
  _polarization(1) *=1e-4;
  _polarization(2) *=1e-4;

  return ( _polarization);

}



//  std::vector< std::map< ID, double > > solution;

//  if  (_strain_sim->get_solution(_elem,q_point,pol_ID,solution))
//   {

//     Tensor2Sym strain;

//     strain(1,1) = solution[0].find(var_map[E_XX])->second;
//     strain(2,1) = solution[0].find(var_map[E_XY])->second;
//     strain(3,1) = solution[0].find(var_map[E_XZ])->second;
//     strain(2,2) = solution[0].find(var_map[E_YY])->second;
//     strain(3,2) = solution[0].find(var_map[E_YZ])->second;
//     strain(3,3) = solution[0].find(var_map[E_ZZ])->second;

//     // _dynamical_matrix =  deformation_potential * strain;

//     //Material* mat = get_material();
//     //const RotatedCrystal&   cr = mat->get_rotated_crystal ();
//     //rotate_to_calculation_system(cr.RotMatrix);


//   }

//  //return  (_pyropolarization + _piezopolarization);


//}
