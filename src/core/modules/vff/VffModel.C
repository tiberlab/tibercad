#include "VffModel.h"

#include "SimulationOptions.h"
#include "Material.h"
#include "Database.h"
#include "Messages.h"
#include "RuntimeException.h"


VffModel*
VffModel::create(const Material* mat, const ModelOptions& options)
{

  std::cout << "SONO STATO CREATO E IL MIO MATERIALE E " << mat->get_name() << std::endl;
  return PhysicalModelInterface::create<VffModel>(_create, _destroy, mat, options);

}
//
//void
//VffModel::prepare_submodels(void)
//{
//
//  ModelOptions opts;
//  opts.set_option("type","constant");
////  create_submodel(_tcm, "thermal_conductivity", opts);
////  create_submodels(_hsm, "heat_source");
//
//}

void VffModel::read_database( )
 {

   const Database& db = get_database();
   db.set_section("elasticity");
   std::cout << "DATABASE IS READ " << std::endl;

   if (get_material()->get_structure() == "wz")
   {
       throw RuntimeException("VffModel: wz structure is not supported yet");
       //     _c13 = db.get("C13", 0.0, true);
       //     _c33 = db.get("C33", 0.0, true);
   }

   _c11 = db.get("C11", 0.0, true);
   _c12 = db.get("C12", 0.0, true);
   _c44 = db.get("C44", 0.0, true);

   db.set_section("lattice");
   _a = db.get("a", 0.0, true);

}

void VffModel::do_init( )
 {
  //NOTE: in database the lattice constant is in nm and the stiffness constants
  //are in GPa, therefore alpha and beta are in N/m without unit conversion
  //The reference distance must be in A, as atomic distances are evaluated in Amstrong
  _alpha = (_c11 + 3.0 * _c12) * (_a  / 4.0);
  _beta = (_c11 - _c12) *  (_a  / 4.0);
  _d = _a * (sqrt(3.0) / 4.0) * 10.0;
  _teta = -0.3333;
  std::cout << "alpha is " << _alpha << std::endl;
 }

