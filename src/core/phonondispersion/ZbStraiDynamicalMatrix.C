// $Id$

#include "ZbStrainDynamicalMatrix.h"
#include "getpot.h"
#include "Material.h"
#include "Database.h"
#include "RotatedCrystal.h"  
#include "PhononModel.h"  
#include "SimulationInterface.h" 

//--------------------------------------------------------//

ZbStrainDynamicalMatrix::ZbStrainDynamicalMatrix(const ModelOptions& options) :
  DynamicalMatrix(options),
  deformation_potential(0)
{
}


void  ZbStrainDynamicalMatrix::read_database(void)
{
 
  const Material* mat = get_material();

  GetPot data((mat->get_database()).get_data_file());

  p_norm = data("p_norm", 0.0);
  q_norm = data("q_norm", 0.0);
  r_norm = data("r_norm", 0.0); 
  w0 = data("w0",0.0);
 
}

//---------------------------------------------------------//



void  ZbStrainDynamicalMatrix::do_init(void)
{

   const ModelOptions& options = get_options();


   p_norm =  options.get_option("p_norm",p_norm);
   q_norm =  options.get_option("q_norm",q_norm);
   r_norm =  options.get_option("r_norm",r_norm);
   w0 = options.get_option("w0",w0);
   w0 = w0 / 8065.6;
   //K11
   double p = p_norm * w0 *w0; 
   deformation_potential(1,1,1,1) = p;
   deformation_potential(2,2,2,2) = p;
   deformation_potential(3,3,3,3) = p;
   //K12
   double q = q_norm * w0 *w0;  
   deformation_potential(2,2,1,1) = q ;
   deformation_potential(3,3,2,2) = q;
   deformation_potential(3,3,1,1) = q;
   //K44
   double r = r_norm * w0 *w0; 
   deformation_potential(2,1,2,1) = r;
   deformation_potential(3,1,3,1) = r;
   deformation_potential(3,2,3,2) = r;

   //read strain simulation

   std::string strain_sim = get_options().get_option("strain_simulation", "");

   _simul = SimulationInterface::find_simulation(strain_sim);

  if ( _simul == NULL)
   throw InitFailedException("Could not find " + strain_sim);

    //Get the variable IDs
    var_map[E_XX]=_simul->get_solution_id("eps_xx");
    var_map[E_XY]=_simul->get_solution_id("eps_xy");
    var_map[E_XZ]=_simul->get_solution_id("eps_xz");
    var_map[E_YY]=_simul->get_solution_id("eps_yy");
    var_map[E_YZ]=_simul->get_solution_id("eps_yz");
    var_map[E_ZZ]=_simul->get_solution_id("eps_zz");

   std::map<ID,ID>::iterator      it(var_map.begin());
   std::map<ID,ID>::iterator      end(var_map.end());
   for(; it!=end; ++it)
     ID_set.insert(it->second);

}

void  ZbStrainDynamicalMatrix::re_init(void)
{


  const Elem* elem = _phonon_model->get_element(); 

  std::vector<Point> h_point(1);
  h_point[0] = elem->vertex_average();

  std::vector< std::map< ID, double > > solution;

  if  (_simul->get_solution(elem,h_point,ID_set,solution))
  {
    
    Tensor2Sym strain;
    
    strain(1,1) = solution[0].find(var_map[E_XX])->second;
    strain(2,1) = solution[0].find(var_map[E_XY])->second;
    strain(3,1) = solution[0].find(var_map[E_XZ])->second;
    strain(2,2) = solution[0].find(var_map[E_YY])->second;
    strain(3,2) = solution[0].find(var_map[E_YZ])->second;
    strain(3,3) = solution[0].find(var_map[E_ZZ])->second;
    
    _dynamical_matrix =  deformation_potential * strain;
    
    Material* mat = get_material();
    
    const RotatedCrystal&   cr = mat->get_rotated_crystal ();
    
    rotate_to_calculation_system(cr.RotMatrix);
    
    
  }
   
    

}
