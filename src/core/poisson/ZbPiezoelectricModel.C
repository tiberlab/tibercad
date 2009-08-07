// $Id$

#include "ZbPiezoelectricModel.h"
#include "Material.h"
#include "Database.h"
#include "RotatedCrystal.h"
#include "elem.h"
#include "tensor.h"
//---------------------------------------------------//

 
//----------------------------------------------//
void ZbPiezoelectricModel::read_database ( )
{
 
  Database& db = get_database();
  db.set_section("piezoelectricity");

  e14 = db.get("e14", 0.0, true);
  
  
}

//-----------------------------------------------//
void ZbPiezoelectricModel::do_init ( )
{
  
  std::string _simul_name = get_options().get_option("strain_simulation", "");
  _simul = SimulationInterface::find_simulation(_simul_name);
  if ( _simul == NULL)
    throw InitFailedException("Could not find " + _simul_name);



  // e14 = get_parameter("e14", e14); 

  var_map.clear();
  var_map[EXX]=_simul->get_variable_id("eps_xx");
  var_map[EXY]=_simul->get_variable_id("eps_xy");
  var_map[EXZ]=_simul->get_variable_id("eps_xz");
  var_map[EYY]=_simul->get_variable_id("eps_yy");
  var_map[EYZ]=_simul->get_variable_id("eps_yz");
  var_map[EZZ]=_simul->get_variable_id("eps_zz");

  std::map<ID,ID>::iterator      it(var_map.begin());
  std::map<ID,ID>::iterator      end(var_map.end());
  for(; it!=end; ++it)
    ID_set.insert(it->second);   

}



void  ZbPiezoelectricModel::calculate_piezopolarization(const Elem* elem)
{


  Material* mat = get_material();

  const RotatedCrystal& cr = mat->get_rotated_crystal ();

  std::vector< std::map< ID, double > > solution;

  std::vector<Point> p(1);

  p[0] = elem->centroid();

   Tensor2Sym eps(0);
  if  (_simul->get_solution(elem,p,ID_set,solution))
  {
    eps(1,1) = solution[0].find(var_map[EXX])->second; 
    eps(2,1) = solution[0].find(var_map[EXY])->second;
    eps(3,1) = solution[0].find(var_map[EXZ])->second;
    eps(2,1) = solution[0].find(var_map[EXY])->second; 
    eps(2,2) = solution[0].find(var_map[EYY])->second;
    eps(3,3) = solution[0].find(var_map[EZZ])->second;


    eps = sym(cr.RotMatrix.transpose() * ( eps * (cr.RotMatrix)));

    _P(1) = 2*e14*eps(3,2);
    _P(2) = 2*e14*eps(3,1);
    _P(3) = 2*e14*eps(2,1);



    rotate_to_calc_system(cr.RotMatrix);
   

  }

  

}


