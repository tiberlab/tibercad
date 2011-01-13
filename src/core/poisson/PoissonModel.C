#include "PoissonModel.h"
#include "Material.h"
#include "ChargeDensityModel.h"
#include "Database.h"
#include "getpot.h"
#include "PiezoelectricModel.h"


PoissonModel::PoissonModel(const ModelOptions& options) :
  PhysicalModel(options),
  _strain_sim(NULL),
  _chd_sim(NULL),
  _epsilon_model(NULL),
  _epsilon(0),
  _piezo_model(NULL)
{
}



PoissonModel::~PoissonModel(void)
{

  destroy(_piezo_model);
  destroy(_epsilon_model);

}

//==========================================================================//

PhysicalModelInterface* PoissonModel::create_new (void) const
{
   return new PoissonModel(get_options());
}

//==========================================================================//
void PoissonModel::do_init()
{
}
  //const ModelOptions& opt =  get_options ();

// <<<<<<< .mine

//    model_opt.piezo_pol = get_options().get_option("piezoelectric_field", false);

//    if  (model_opt.piezo_pol == true)
//    {
//      std::string strain_simul = get_options().get_option("strain_simulation", "no");
//      _strain_if.set_simulation(strain_simul);
//    }

//   //Pyropolarization
//   model_opt.pyro_pol = get_options().get_option("piroelecric_field", false);
//   if (model_opt.pyro_pol)
//   {
//     PhysicalModelInterface::destroy(_pyropolarization);
//     _pyropolarization = PyroPolarization::create(get_material());
//     _pyropolarization->set_material(get_material());
//     _pyropolarization->set_simulator_id(get_simulator_id());
//     _pyropolarization->init();
//   }

// =======
// >>>>>>> .r1571

   //Density Charge Model

  //    model_opt.chd_sim_name = get_options().get_option("charge_density_simulation", "no_sim");

      //if  (model_opt.chd_sim_name.compare("no_sim")==1)
      //{
//	_chd_sim =  SimulationInterface::find_simulation( model_opt.chd_sim_name);

//	if (_chd_sim == NULL)
//	  throw InitFailedException("Unknown charge_density simulation: " +   model_opt.chd_sim_name );
//	else
//	  charge_id = _chd_sim->get_variable_id("charge_density");
  //    }
      //doping

    //  model_opt.add_doping = get_options().get_option("add_doping",false);
      //---------------



void PoissonModel::create_submodels()
{

  ModelOptions::const_submodel_iterator it;
  ModelOptions::const_submodel_iterator end;


  //Piezopolarization-----
  destroy(_piezo_model);
  it = get_options().submodels_begin("piezo_polarization");
  end = get_options().submodels_end("piezo_polarization");

  if (it != end)
  {

    _piezo_model = dynamic_cast<PiezoelectricModel*>(
                     PhysicalModelInterface::create("piezoelectric_model_" +
                    get_material()->get_structure(), it->second));

   if (_piezo_model == NULL)
      throw InitFailedException("Could not create piezoelectric model");

   _piezo_model->set_owner(get_material());
   _piezo_model->init();

  }



   //Dielectric constant
   destroy(_epsilon_model);

   const ModelOptions& opt =  get_options ();
   _epsilon_model = OptDielectricConstant::create(get_material()->get_structure(), opt);
   _epsilon_model->set_owner(get_material());
   _epsilon_model->init();
   _epsilon_model->get_dielectric_real(_epsilon);


 }



//==========================================================================//
void PoissonModel::do_init_alloy (const PhysicalModelInterface *comp_A, const PhysicalModelInterface *comp_B, double xa)
{

  const PoissonModel* matA = dynamic_cast< const PoissonModel*> (comp_A);
  const PoissonModel* matB = dynamic_cast< const PoissonModel*> (comp_B);

  destroy(_epsilon_model);
 // dielectric_model = create_submodel_alloy(matA->dielectric_model,matB->dielectric_model,xa);

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


void
PoissonModel::get_dielectric_constant(RealTensor& epsilon)
{


  for (unsigned int i=0; i<3;i++)
    for (unsigned int j=0; j<=i;j++)
    {
      epsilon(i,j) = _epsilon(i+1,j+1);
      epsilon(j,i) = epsilon(i,j);
    }



  epsilon *=  Constants::epsilon * 1e-2;

}


void
PoissonModel::get_total_polarization(std::vector<RealGradient>& p,const std::vector< Point >& points)
{

  std::vector<RealGradient> piezo_pol(points.size());
  get_piezo_polarization(piezo_pol,points);

  std::vector<RealGradient> pyro_pol(points.size());
  get_pyro_polarization(pyro_pol,points);

  for (ID n = 0; n < points.size();n ++)
    p[n] = piezo_pol[n] + pyro_pol[n];

}

void
PoissonModel::get_piezo_polarization(std::vector<RealGradient>& p,const std::vector< Point >& points)
{
 std::vector<RealGradient> _polarization(points.size());

  if (_piezo_model != NULL)
  {

    for (ID n = 0; n< points.size(); n++)
      {
	const Point p = points[n];
	_piezo_model->calculate_piezopolarization(_elem,p);

	Tensor1 piezo(0);
	_piezo_model->get_piezopolarization(piezo);

        _polarization[n](0);
	_polarization[n](0) = piezo(1);
	_polarization[n](1) = piezo(2);
	_polarization[n](2) = piezo(3);
      }

  }

  for (ID n = 0; n< points.size();n++)
    {
      //Scaling from Q/m2 tp Q/cm2
      _polarization[n](0) *=1e-4;
      _polarization[n](1) *=1e-4;
      _polarization[n](2) *=1e-4;
    }

  p = _polarization;


}

void
PoissonModel::get_pyro_polarization(std::vector<RealGradient>& p,const std::vector< Point >& points)
{

 std::vector<RealGradient> _polarization(points.size());


  for (ID n = 0; n< points.size();n++)
    {
      //Scaling from Q/m2 tp Q/cm2
      _polarization[n](0) *=1e-4;
      _polarization[n](1) *=1e-4;
      _polarization[n](2) *=1e-4;
    }

  p = _polarization;

}
