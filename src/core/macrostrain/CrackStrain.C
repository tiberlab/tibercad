#include "CrackStrain.h"
#include "SimulationEnvironment.h"
#include "mesh.h"
#include "MacrostrainModel.h"
#include "Material.h"

using namespace std;



void CrackStrain::parse_options()
{
  

  const ModelOptions& opt = get_options();

  _x0 = opt.get_option("x0", 0.0);

  _y0 = opt.get_option("y0", 0.0);

  // units: GPa m^1/2
  _Ki = opt.get_option("Ki", 0.0);

  _sigma_ys = opt.get_option("sigma_ys", 0.0);

 


}

//------------------------------------------------------//

void CrackStrain::do_init()
{

  StrainSimulation::do_init();

}

//------------------------------------------------------//


void CrackStrain::do_solve()
{

  parse_options();

  result_strain.clear();

  const MeshBase& mesh = get_environment().get_mesh();

  Tensor2Sym stress;

  Tensor2Sym strain;
  

  MeshBase::const_element_iterator       el     = mesh.active_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_elements_end();

  for ( ; el != end_el ; ++el) 
  {
    const Elem* elem = *el;
    const Point p = elem->centroid();

    calculate_stress(stress, p(0), p(1)); 

    ID subdomain = elem->subdomain_id();
    
    const Material* mat = _device->get_material(subdomain);

    MacrostrainModel* macrostrain_model = dynamic_cast<MacrostrainModel*>(   mat ->get_model(get_id())     );

    const Stiffness* stiffness = macrostrain_model->get_stiffness();

    stiffness->calculate_strain_from_stress (stress,  strain);



    const RotatedCrystal* crystal_el = &(mat->get_rotated_crystal());
    
    Tensor2Gen RotM = (crystal_el->RotMatrix).transpose();//get rotation matrix

    // RotM  = 1;
      
    Tensor2Gen eps1 = (RotM * strain) * RotM.transpose();  //transform to crystal system

    strain  = sym(eps1); //result has to be symmetric

    assert (::norm(strain - eps1) < 1e-6); //is it really symmetric


    result_strain.insert(std::pair <const Elem*, Tensor2Sym>  (elem, strain));

    
    // result_strain.insert(std::pair <const Elem*, Tensor2Sym>  (elem, stress));
    
  }

  
}

//-------------------------------------------------------//

void CrackStrain::calculate_stress(Tensor2Sym& stress, const double x_input, const double y_input) const
{

  stress = 0;

  double mesh_units =  _device->get_mesh_units();

 

  double x = (x_input - _x0)*mesh_units;
  double y = (y_input - _y0)*mesh_units;

  double rho = std::sqrt(x*x + y*y);
  double teta;

  if (rho < 1e-14)
  {
    rho = 0.0;
    teta = 0.0;
  }
  else
  {
   
    if (x>0) 
      teta = std::atan2(y,x);
    else if ((x<0)&&(y>=0))
      teta = std::atan2(y,x) + M_PI;
    else if ((x<0) && (y < 0))
      teta = std::atan2(y,x) - M_PI;
    else if (x == 0 && y > 0)
      teta = M_PI/2.0;
    else 
      teta = -M_PI/2.0;
  }



  //rho_b *= std::pow(std::cos(teta/2.0) * (1 + std::abs(std::sin(teta/2.0))), 2.0);
  // Tresca
  //double rho_b = _Ki * _Ki / (2 * M_PI * _sigma_ys * _sigma_ys);
  //rho_b *= std::pow(std::cos(teta/2.0) * (1 + std::sin(std::abs(teta)/2.0)), 2.0);
  // Von Mises
  double rho_b = _Ki * _Ki / (4 * M_PI * _sigma_ys * _sigma_ys);
  rho_b *= 1 + 1.5 * std::sin(teta) * std::sin(teta) + std::cos(teta);


  //if (rho < rho_b)
  if (rho < 2 * rho_b)
  {
    stress(1,1) = _sigma_ys;
    stress(2,2) = _sigma_ys;
    stress(2,1) = 0;

    if (x < 0.0)
    {
      stress(1,1) = -_sigma_ys;
      stress(2,2) = -_sigma_ys;
    }
  }
  else
  {

    //double t = _Ki/std::sqrt(2.0*M_PI*(rho));
    double t = _Ki/std::sqrt(2.0*M_PI*(rho - rho_b));

    stress(1,1) = t * std::cos(teta/2.0) * (1.0 - std::sin(teta/2.0)*std::sin(3.0/2.0*teta));
 
    stress(2,2) = t * std::cos(teta/2.0) * (1.0 + std::sin(teta/2.0)*std::sin(3.0/2.0*teta));
  
    stress(2,1) = t * std::sin(teta/2.0) * std::cos(teta/2.0) * std::cos(3.0*teta/2.0);
  }


}

//----------------------------------------------------------------------------------------------//

void CrackStrain::build_elemental_results(const std::set<std::string>& variables,
				       std::vector<double>& results, std::vector<std::string>& legend) 
{

  std::vector<std::string> eps_names;
  std::vector<double> eps_data;  

  std::vector<std::string> stress_names;
  std::vector<double> stress_data;  

  std::vector<std::string> pol_names;
  std::vector<double> pol_data;  

  prepare_strain_data_for_output( eps_names,  eps_data);
  prepare_stress_data_for_output( stress_names,  stress_data);
  prepare_polarization_data_for_output( pol_names,  pol_data);

  // for (unsigned int i= 0; i < pol_data.size() ; i++)
  //  cerr << pol_data[i] << "\n";


  
  const set<string>::const_iterator varend = variables.end();
  
  const string strain_name("strain");
  const string stress_name("stress");
  const string pol_name("polarization");

  int num_var = 0;

  int strain_id = -1;
  if (variables.find(strain_name) != varend)
  {
    strain_id = num_var;
    num_var += 6;
  }

  int stress_id = -1;
  if (variables.find(stress_name) != varend)
  {
    stress_id = num_var;
    num_var += 6;  
  }

  int pol_id = -1;
  if (variables.find(pol_name) != varend)
  {
    pol_id = num_var;
    num_var += 3;
  }

  unsigned int num_elem = eps_data.size()/6;



  results.resize(num_var * num_elem);
  legend.resize(num_var);

  

  if (strain_id != -1)
  {//we do strain
    for (short i = 0; i < 6; i++)
    {	
      legend[i] = eps_names[i];
      for (unsigned int j = 0; j < num_elem; j++)
	results[i + j * num_var ] = eps_data[i + j * 6];  
    }
    
  }
   

  if (stress_id != -1)
  {//we do stress
    for (short i = 0; i < 6; i++)
    {	
      legend[i + stress_id] = stress_names[i];
      for (unsigned int j = 0; j < num_elem; j++)
	results[i + j * num_var + stress_id] = stress_data[i + j * 6];  
    }
    
  }
  


  
  if (pol_id != -1)
  {//now we do polarization
    for (short i1 = 0; i1 < 3; i1++)
    {  
      
      legend[i1 + pol_id] = pol_names[i1];
      for (unsigned int j = 0; j < num_elem; j++)
	results[i1 + j * num_var + pol_id] = pol_data[i1 + j * 3];  
      
      
    }
  }



}


//---------------------------------------------------------------------------//
void CrackStrain::prepare_strain_data_for_output( std::vector<std::string>& eps_names, std::vector<double>& eps_data )
{

  const MeshBase& mesh = get_environment().get_mesh();

  char num_i[2];
  char num_j[2];
  string eps_ij;

  
  unsigned int Number_of_elements = mesh.n_active_elem();

  eps_data.resize(Number_of_elements*6);
  eps_names.resize(6);


  unsigned int index = 0;
  for (int i = 1; i <=3 ; i++)
    for (int j = 1; j <=i; j++)
    {

      sprintf( num_i, "%i",i);
      sprintf( num_j, "%i",j);
	

      eps_ij = "eps_" + string(num_i) + string(num_j);

      eps_names[index] = eps_ij ;

      index++;

    }

  


  


  MeshBase::const_element_iterator       el     = mesh.active_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_elements_end();

  unsigned int elem_number = 0;
      
  Tensor2Sym strain;
      
  for ( ; el != end_el ; ++el) 
  {

    const Elem* elem = *el;

    ID subdomain = elem->subdomain_id();
    
    const Material* mat = _device->get_material(subdomain);

	
    const RotatedCrystal* crystal_el = &(mat->get_rotated_crystal());

    Tensor2Gen RotM = crystal_el->RotMatrix;//get rotation matrix
    
    Tensor2Gen eps1 = (RotM * result_strain[elem]) * RotM.transpose();  //transform to calculation system
    
    strain  = sym(eps1); //result has to be symmetric

    index = 0;
	
    for (int i = 1; i <=3 ; i++)
      for (int j = 1; j <=i; j++)
      {

	eps_data[index + elem_number * 6  ] = strain(i,j);
	    
	index++;
      }


    elem_number++;
  }
  
      


}



//--------------------------------------------------------------------------//

void CrackStrain::prepare_polarization_data_for_output( std::vector<std::string>& polariz_names, std::vector<double>& polariz_data )
{
  char num_i[2];

  const MeshBase& mesh = get_environment().get_mesh();


  unsigned int Number_of_elements = mesh.n_active_elem();

  polariz_data.resize(Number_of_elements*3);

  polariz_names.resize(3);

  polariz_names[0] = "Px";
  polariz_names[1] = "Py";
  polariz_names[2] = "Pz";


  MeshBase::const_element_iterator       el     = mesh.active_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_elements_end();

  unsigned int el_number = 0;

  Tensor1 polariz_vec;

  for ( ; el != end_el ; ++el) 
  {
      
    Elem* elem = *el;

    ID subdomain = elem->subdomain_id();

    const Material* mat = _device->get_material(subdomain);

    const RotatedCrystal* crystal_el = &(mat->get_rotated_crystal());

    MacrostrainModel* macrostrain_model =  dynamic_cast<MacrostrainModel*>(   mat ->get_model(get_id())     );

    Tensor1 polariz = ( macrostrain_model->get_piezo() ) -> get_polariz_cryst(result_strain[elem]); //crystal system

   


    polariz_vec =(crystal_el->RotMatrix) * polariz; //calculation system

   

   
    
    polariz_data[0 + el_number*3] = polariz_vec (1);
    polariz_data[1 + el_number*3] = polariz_vec (2);
    polariz_data[2 + el_number*3] = polariz_vec (3);

    el_number++;

  }



}




void
CrackStrain::prepare_stress_data_for_output( std::vector<std::string>& stress_names, std::vector<double>& stress_data ) 
{

  const MeshBase& mesh = get_environment().get_mesh();

  char num_i[2];
  char num_j[2];
  string stress_ij;

  
  unsigned int Number_of_elements = mesh.n_active_elem();

  stress_data.resize(Number_of_elements*6);
  stress_names.resize(6);


  unsigned int index = 0;
  for (int i = 1; i <=3 ; i++)
    for (int j = 1; j <=i; j++)
    {

      sprintf( num_i, "%i",i);
      sprintf( num_j, "%i",j);
	

      stress_ij = "stress_" + string(num_i) + string(num_j);

      stress_names[index] = stress_ij ;

      index++;

    }

  


  


  MeshBase::const_element_iterator       el     = mesh.active_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_elements_end();

  unsigned int elem_number = 0;
      
  Tensor2Sym strain;
      
  for ( ; el != end_el ; ++el) 
  {

    const Elem* elem = *el;

    ID subdomain = elem->subdomain_id();
    
    const Material* mat = _device->get_material(subdomain);
    MacrostrainModel* macrostrain_model =
      dynamic_cast<MacrostrainModel*>(mat->get_model(get_id())); 

	
    const RotatedCrystal* crystal_el = &(mat->get_rotated_crystal());

    Tensor2Gen RotM = crystal_el->RotMatrix;//get rotation matrix
    
    Tensor2Gen eps1 = (RotM * result_strain[elem]) * RotM.transpose();  //transform to calculation system
    
    strain  = sym(eps1); //result has to be symmetric

    //Elasticity in the calculation system
    Stiffness* C_tensor_el;
    C_tensor_el = macrostrain_model->get_stiffness();
    Tensor4DSym C_calc =  C_tensor_el->C_calc;

    Tensor2Sym stress_el = strain * C_calc;       


    index = 0;
	
    for (int i = 1; i <=3 ; i++)
      for (int j = 1; j <=i; j++)
      {

	stress_data[index + elem_number * 6  ] = stress_el(i,j);
	    
	index++;
      }


    elem_number++;
  }
}
