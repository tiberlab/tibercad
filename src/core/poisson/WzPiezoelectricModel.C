// $Id$

#include "WzPiezoelectricModel.h"
#include "Material.h"
#include "Database.h"
#include "RotatedCrystal.h"
#include "elem.h"
#include "tensor.h"
//---------------------------------------------------//


//----------------------------------------------//
void WzPiezoelectricModel::read_database ( )
{

  Database& db = get_database();
  db.set_section("piezoelectricity");

  e33 = db.get("e33", 0.0, true);
  e31 = db.get("e31", 0.0, true);
  e15 = db.get("e15", 0.0, true);



}

//-----------------------------------------------//
void WzPiezoelectricModel::do_init ( )
{

  std::string _simul_name = get_options().get_option("strain_simulation", "");
  _simul = SimulationInterface::find_simulation(_simul_name);
  if ( _simul == NULL)
    throw InitFailedException("Could not find " + _simul_name);



  // e33 = get_parameter("e33", e33);
//   e31 = get_parameter("e31", e31);
//   e15 = get_parameter("e15", e15);


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



void  WzPiezoelectricModel::calculate_piezopolarization(const Elem* elem)
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
    eps(2,2) = solution[0].find(var_map[EYY])->second;
    eps(3,2) = solution[0].find(var_map[EYZ])->second;
    eps(3,3) = solution[0].find(var_map[EZZ])->second;

    //calculation system -> crystal system
    eps = sym(cr.RotMatrix.transpose() * ( eps * (cr.RotMatrix)));


    _pol(1) = 2.0 * e15 * eps(3,1);
    _pol(2) = 2.0 * e15 * eps(3,2);
    _pol(3) = e31 * eps(1,1) + e31*eps(2,2) + e33 * eps(3,3);

    //calculation system -> crystal system
    rotate_to_calc_system(cr.RotMatrix);


  }



}


