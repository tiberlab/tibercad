#include "PhononModel.h"
#include "Material.h"
#include "RamanTensor.h"
void PhononModel::fake(void)
{
}


PhononModel::PhononModel() :
  _elem(NULL)
{
}

PhononModel::~PhononModel()
{
   clear_dynamical_matrix_models();

   PhysicalModelInterface::destroy(raman_tensor_model);
   
}

void PhononModel::get_full_dynamical_matrix(Tensor2Sym& dynamical_matrix)
{

  dyn_mat_iterator it  =  _dynamical_matrix_models.begin();
  dyn_mat_iterator it_end =  _dynamical_matrix_models.end();
  
  
  Tensor2Sym total(0);
  for ( ; it != it_end; ++it)
  {
    Tensor2Sym partial(0);

    (it->second)->get_dynamical_matrix(partial);

    total = total + partial;
    
  }

  dynamical_matrix = total;

}

	
//void PhononModel::get_free_dynamical_matrix(std::vector<std::vector< double > >& D)
void PhononModel::get_free_dynamical_matrix(Tensor2Sym& dynamical_matrix)
{

  //Tensor2Sym DM;
   _dynamical_matrix_models[free_ID]->get_dynamical_matrix(dynamical_matrix);
   //(_dynamical_matrix_models.begin()->second)->get_dynamical_matrix(dynamical_matrix);

}

void PhononModel::get_raman_tensor(std::vector<Tensor2Sym>& raman_tensor)
{

 raman_tensor_model->get_raman_tensor(raman_tensor);

}


void PhononModel::do_init(void)
{

  const ModelOptions dummy;
  free_ID = add_dynamical_matrix_model("free_dynamical_matrix_" +
			     get_material()->get_structure(),dummy);

  ModelOptions::const_submodel_iterator it;
  ModelOptions::const_submodel_iterator end;
  
  //Heat source models
  it = get_options().submodels_begin("dynamical_matrix");
  end = get_options().submodels_end("dynamical_matrix");

  for ( ; it != end; ++it)
  {
    const std::string& name = (it->second).get_option("model", "");
  
    ID id = add_dynamical_matrix_model(name + "_" +
			       get_material()->get_structure(),it->second);
  }


  std::string model_name =  "raman_tensor_" +get_material()->get_structure();      

  raman_tensor_model = dynamic_cast<RamanTensor*>(PhysicalModelInterface::create(model_name)); 
  
  if (raman_tensor_model == NULL)
    throw InitFailedException("No such raman tensor model" + model_name);
  
  raman_tensor_model->set_phonon_model(this);
  raman_tensor_model->set_material(get_material());
  raman_tensor_model->set_simulator_id(get_simulator_id());
  raman_tensor_model->init();


  const ModelOptions opt = get_options();

  std::vector<double> LP(3,0.0);
  LP[0] = 1.0;
  LP[1] = 0.0;
  LP[2] = 0.0;
  opt.get_option("LightPolarization",LP);

  std::vector<double> LD(3,0.0);
  LD[0] = 0.0;
  LD[1] = 0.0;
  LD[2] =-1.0;
  opt.get_option("LightDirection",LD);
  


    //vectorProduct()
    _light_polarization.resize(4);

    _light_polarization[0](1) = LP[0] ;   
    _light_polarization[0](2) = LP[1] ;
    _light_polarization[0](3) = LP[2] ;

    _light_polarization[1](1) = LP[0] ;   
    _light_polarization[1](2) = LP[1] ;
    _light_polarization[1](3) = LP[2] ;

     Tensor1  dten;   
     dten(1) = LD[0];
     dten(2) = LD[1];
     dten(3) = LD[2];

    _light_polarization[2] = vectorProduct(_light_polarization[0],dten);
 
    _light_polarization[3] =  _light_polarization[2] +   _light_polarization[1];   


  
}


ID
PhononModel::add_dynamical_matrix_model(const std::string& model_name, 
                                const ModelOptions& options)
{
  
  DynamicalMatrix* model = dynamic_cast<DynamicalMatrix*>(
							  PhysicalModelInterface::create(model_name,options)); 
  
  if (model == NULL)
    throw InitFailedException("No such dynamical matrix model" + model_name);
  
  model->set_phonon_model(this);
  model->set_material(get_material());
  model->set_simulator_id(get_simulator_id());
  model->init();
  
  ID id = model->get_id();
  _dynamical_matrix_models[id] = model;

  
  return id;   

}
void
PhononModel::clear_dynamical_matrix_models(void)
{

dyn_mat_iterator it =   _dynamical_matrix_models.begin();
dyn_mat_iterator end =  _dynamical_matrix_models.end();

for ( ; it != end; ++it)
 {
     PhysicalModelInterface::destroy(it->second);
 }

 _dynamical_matrix_models.clear();
    

}

void
PhononModel::get_light_polarization(std::vector<Tensor1>& light_polarization)
{

  light_polarization = _light_polarization;

}

//void
//PhononModel::re_init(void)
//{
//  update_dynamical_matrix();
//}





  //  std::cerr<<DM<<std::endl;

  //unsigned int n = 3;   
  //D.resize(n);
  //for (unsigned int i = 0; i<n; i++)
  //{
  // D[i].resize(n,0.0);
  // for (unsigned int j = 0; j<n; j++)
  // {
  //   double value;
  //   if (i < j) 
  //	value = DM(j+1, i+1);
  //    else
  //	value = DM(i+1, j+1);	   
      
  //    D[i][j] = value;
  //  }
  //} 
  
  //D[0][0] = 1.0;
  // D[0][1] = 3.0;
  // D[1][0] = 3.0;
  // D[1][1] = 1.0; 

//      D = _dynamical_matrix;



//PhononModel::~PhononModel()
//{
//}

//==========================================================================//

// PhysicalModelInterface* PhononModel::create_new (void) const
// {
//   return new PhononModel();
// }

// //==========================================================================//
// void PhononModel::do_init()
// {

//   ModelOptions::const_submodel_iterator it;
//   ModelOptions::const_submodel_iterator end;
  
//   //Heat source models
//   it = get_options().submodels_begin("dynamical_tensor");
//   end = get_options().submodels_end("dynamical_tensor");

//   for ( ; it != end; ++it)
//   {
//     //const std::string& name = (it->second).get_option("model", "");
    
//   }



//   if (it != end)
//   {
    
    
//     //kappa =dynamic_cast<LatticeThermalConductivity*>(
//     //	      PhysicalModelInterface::create("lat_therm_cond_" +
//     //	      get_material()->get_structure(), it->second)); 
    
  
//     //if (kappa == NULL)
//     // throw InitFailedException("Could not create lattice thermal conductivity model");

 
//     }
//     else
//     {
//     //kappa = dynamic_cast<LatticeThermalConductivity*>(
//     //  PhysicalModelInterface::create("lat_therm_cond_" +
//     //	get_material()->get_structure()));

//      }
//     //kappa->set_temperature(SimulationOptions::temperature);
    
//     //kappa->set_material(get_material());

//     //kappa->init();




   
// }



// //==========================================================================//
// void PhononModel::do_init_alloy (const PhysicalModelInterface *comp_A, const PhysicalModelInterface *comp_B, double xa)
// {

//   const PhononModel* matA = dynamic_cast< const PhononModel*> (comp_A);

//   const PhononModel* matB = dynamic_cast< const PhononModel*> (comp_B);

//   //kappa->init_alloy(matA->kappa, matB->kappa, xa);

  

// }
