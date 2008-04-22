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

  _Ki = opt.get_option("Ki", 0.0);

 


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

  const Mesh& mesh = get_environment().get_mesh();

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
      
    Tensor2Gen eps1 = (RotM * strain) * RotM.transpose();  //transform to crystal system

    strain  = sym(eps1); //result has to be symmetric

    assert (norm(strain - eps1) < 1e-6); //is it really symmetric

    result_strain.insert(std::pair <const Elem*, Tensor2Sym>  (elem, strain));

    
    //result_strain.insert(std::pair <const Elem*, Tensor2Sym>  (elem, stress));
    
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
    if (x >= 0.0)
      teta = std::asin(y/rho);
    else
      teta = M_PI - std::asin(y/rho);
  }


  double t = _Ki/std::sqrt(2*M_PI*rho);

 

  stress(1,1) = t * std::cos(teta/2) * (1 - std::sin(teta/2)*std::sin(3/2*teta));
 
  stress(2,2) = t * std::cos(teta/2) * (1 + std::sin(teta/2)*std::sin(3/2*teta));
  
  stress(2,1) = t * std::sin(teta/2) * std::cos(teta/2) * std::cos(3*teta/2);
  

}

//----------------------------------------------------------------------------------------------//

void CrackStrain::build_elemental_results(const std::set<std::string>& variables,
				       std::vector<double>& results, std::vector<std::string>& legend) 
{

  std::vector<std::string> eps_names;
  std::vector<double> eps_data;  

  std::vector<std::string> pol_names;
  std::vector<double> pol_data;  

   prepare_strain_data_for_output( eps_names,  eps_data);
   prepare_polarization_data_for_output( pol_names,  pol_data);


  
  short num_var = 0;
  const set<string>::const_iterator varend = variables.end();
  
  const string strain_name("strain");
  const string pol_name("polarization");

  if (variables.find(strain_name) != varend) num_var += 6;  
  if (variables.find(pol_name) != varend) num_var +=3;

  unsigned int num_elem = eps_data.size()/6;



  results.resize(num_var * num_elem);
  legend.resize(num_var);

  

  if (variables.find(strain_name) != varend)
  {//we do strain
    for (short i = 0; i < 6; i++)
    {	
      legend[i] = eps_names[i];
      for (unsigned int j = 0; j < num_elem; j++)
	results[i + j * num_var ] = eps_data[i + j * 6];  
    }
    
  }
  

  
  if (variables.find(pol_name) != varend)
  {//now we do polarization
    for (short i1 = 0; i1 < 3; i1++)
    {  
      
      legend[i1 + num_var - 3] = pol_names[i1];
      for (unsigned int j = 0; j < num_elem; j++)
	results[i1 + num_var - 3 + j * num_var ] = pol_data[i1 + j * 3];  
      
      
    }
  }



}


//---------------------------------------------------------------------------//
void CrackStrain::prepare_strain_data_for_output( std::vector<std::string>& eps_names, std::vector<double>& eps_data )
{

  const Mesh& mesh = get_environment().get_mesh();

  char num_i[2];
  char num_j[2];
  string eps_ij;

  unsigned int index = 0;
  unsigned int Number_of_elements = mesh.n_active_elem();

  eps_data.resize(Number_of_elements*6);
  eps_names.resize(6);

  for (int i = 1; i <=3 ; i++)
    for (int j = 1; j <=i; j++)
    {

      sprintf( num_i, "%i",i);
      sprintf( num_j, "%i",j);
	

      eps_ij = "eps_" + string(num_i) + string(num_j);

      eps_names[index] = eps_ij ;


      MeshBase::const_element_iterator       el     = mesh.active_elements_begin();
      const MeshBase::const_element_iterator end_el = mesh.active_elements_end();

      unsigned int elem_number = 0;

      for ( ; el != end_el ; ++el) 
      {

	const Elem* elem = *el;

	
	

	eps_data[index + elem_number * 6  ] = result_strain[elem] (i,j);
	elem_number++;
      }
  
      index++;

    }
}



//--------------------------------------------------------------------------//

void CrackStrain::prepare_polarization_data_for_output( std::vector<std::string>& polariz_names, std::vector<double>& polariz_data )
{
  char num_i[2];

  const Mesh& mesh = get_environment().get_mesh();


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

    polariz =(crystal_el->RotMatrix) * polariz; //calculation system
    
    polariz_data[0 + el_number*3] = polariz_vec (1);
    polariz_data[1 + el_number*3] = polariz_vec (2);
    polariz_data[2 + el_number*3] = polariz_vec (3);

    el_number++;

  }



}
