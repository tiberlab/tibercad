#include "VffModel.h"

#include "SimulationOptions.h"
#include "Material.h"
#include "Database.h"
#include "Messages.h"
#include "RuntimeException.h"




VffModel*
VffModel::create(const Material* mat, const ModelOptions& options)
{

  return PhysicalModel::create<VffModel>(_create, _destroy, mat, options);

}
//
void
VffModel::prepare_submodels(void)
{

  ModelOptions opts;
  opts.set_option("type", "automatic");
  create_submodel(_keating, "keating", opts);

}

void
VffModel::read_database( )
{

}

void
VffModel::do_init( )
{
  //NOTE: in database the lattice constant is in nm and the stiffness constants
  //are in GPa, therefore alpha and beta are in N/m without unit conversion
  //The reference distance must be in A, as atomic distances are evaluated in Amstrong
}


double
VffModel::get_alpha(const Atom& atm1, const Atom& atm2) const
{
  if (get_material()->get_structure() == "zb")
    return _keating->get_alpha_0();
  else if (along_c(atm1, atm2))
    return _keating->get_alpha_1();
  else
    return _keating->get_alpha_0();
}


double
VffModel::get_beta(const Atom& atm1, const Atom& atm2, const Atom& atm3) const
{
  if (get_material()->get_structure() == "zb")
    return _keating->get_beta_0();
  else if (along_c(atm1, atm2)  || along_c(atm1, atm3) || along_c(atm2, atm3))
    return _keating->get_beta_1();
  else
    return _keating->get_beta_0();
}


double
VffModel::get_d(const Atom& atm1, const Atom& atm2) const
{
  if (get_material()->get_structure() == "zb")
    return _keating->get_d_0();
  else if (along_c(atm1, atm2))
    return _keating->get_d_1();
  else
    return _keating->get_d_0();
}


double
VffModel::get_costeta(const Atom& atm1, const Atom& atm2, const Atom& atm3) const
{
  if (get_material()->get_structure() == "zb")
    return _keating->get_costeta_0();
  else if (along_c(atm1, atm2) || along_c(atm1, atm3) || along_c(atm2, atm3))
  {
    return _keating->get_costeta_1();
  }
  else
  {
    return _keating->get_costeta_0();
  }
}


bool
VffModel::along_c(const Atom& atm1, const Atom& atm2) const
{
  //Note: tolerance must ne high because it's supposed to work even when the
  //material is strained
  double tol = 0.2;
  double x_d = atm1.get_position(0) - atm2.get_position(0);
  double y_d = atm1.get_position(1) - atm2.get_position(1);
  double z_d = atm1.get_position(2) - atm2.get_position(2);


  //TODO: doing this operation on the fly is not efficient at all, as the
  //c-axis has the same direction in the whole material. You could move it in
  //Material. The loss is not dramatic as this operation is done only once during the
  //parameters matrix assembly, and the operation itself is quite fast
  const libMesh::RealTensor& RotM = get_material()->get_rotation_matrix();
  libMesh::RealVectorValue zz(0); zz(2) = 1.0;
  libMesh::RealVectorValue cc = RotM * zz;
  libMesh::RealVectorValue bond_direction(0);
  bond_direction(0) = x_d; bond_direction(1) = y_d; bond_direction(2) = z_d;
  double scalar_product = pow((cc * bond_direction), 2);
  double norm_product = (cc * cc) * (bond_direction * bond_direction);
  double diff = scalar_product - norm_product;

  return(fabs(diff) < tol);
}


void
VffModel::do_print_info(void)
{
  _keating->print_info();
}
