#include "StrainSimulation.h"
#include "MacrostrainModelInterface.h"
#include "ModelOptions.h"
#include "SimulationEnvironment.h"
#include "MacrostrainModel.h"
#include "Material.h"

Device*  StrainSimulation:: _device;

PhysicalModel*
StrainSimulation::create_physical_model(const ModelOptions& options,
    const Material* mat) const throw (ModelErrorException)
{
  const std::string& modelname = options.get_option("model", "macrostrain");

  MacrostrainModelInterface* model =
    MacrostrainModelInterface::create(modelname, options);

  if (model == NULL)
    throw ModelErrorException(
        "StrainSimulation: No such physical model: " + modelname);

  return model;
}


ID StrainSimulation::convert_variable_name_to_id(const std::string& variable_name) const
{
  
  ID id = INVALID_ID;
 
  if (variable_name == "") return id;

  if (variable_name == "eps_xx")
    id = EPS_XX;
  else if  (variable_name == "eps_yy")
    id = EPS_YY;
  else if  (variable_name == "eps_zz")
    id = EPS_ZZ;
  else if  (variable_name == "eps_xy" || variable_name == "eps_yx")
    id = EPS_XY;
  else if  (variable_name == "eps_yz" || variable_name == "eps_zy")
    id = EPS_YZ;
  else if  (variable_name == "eps_xz" || variable_name == "eps_zx")
    id = EPS_XZ;
  else if  (variable_name == "Px")
    id = P_X;
  else if  (variable_name == "Py")
    id = P_Y;
  else if  (variable_name == "Pz")
    id = P_Z;

  return id;
  
}


void StrainSimulation::get_solution_secure(const Elem* elem,
				   const std::vector<Point>& p, const std::set<ID>& ids,
				   std::vector<std::map<ID, double> >& values)
{
  unsigned int np = p.size();
 
 

  Tensor2Sym strain_el = result_strain[elem];

  ID subdomain = elem->subdomain_id();
          
  const Material* mat = _device->get_material(subdomain);
  
  const RotatedCrystal* crystal_el = &(mat->get_rotated_crystal());

  MacrostrainModel* macrostrain_model =  dynamic_cast<MacrostrainModel*>(   mat ->get_model(get_id())     );
      
  Tensor1 polariz = ( macrostrain_model->get_piezo() )-> get_polariz_cryst(strain_el); //crystal system

  polariz = (crystal_el->RotMatrix) * polariz; //calculation system
 
  for (unsigned int n = 0; n < np; n++)
  {
    if (ids.count(EPS_XX))
      values[n][EPS_XX] = strain_el(1,1);

    if (ids.count(EPS_YY))
      values[n][EPS_YY] = strain_el(2,2);

    if (ids.count(EPS_ZZ))
      values[n][EPS_ZZ] = strain_el(3,3);


    if (ids.count(EPS_XY))
      values[n][EPS_XY] = strain_el(2,1);

    if (ids.count(EPS_YZ))
      values[n][EPS_YZ] = strain_el(3,2);

    if (ids.count(EPS_XZ))
      values[n][EPS_XZ] = strain_el(3,1);
    

    if (ids.count(P_X))
      values[n][P_X] = polariz(1);


    if (ids.count(P_Y))
      values[n][P_Y] = polariz(2);

    if (ids.count(P_Z))
      values[n][P_Z] = polariz(3);


  }
  

}


void StrainSimulation::do_init(void)
{
  SimulationEnvironment& si = get_environment();   

  _device = &( si.get_device() );
}
