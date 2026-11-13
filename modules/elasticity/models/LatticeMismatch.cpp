/*  
 * This file is part of the tiberCAD module elasticity.
 *
 * tiberCAD modules are licensed under the GNU General Public License v3.
 *
 * tiberCAD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published
 * by the Free Software Foundation, either version 3 of the License,
 * or (at your option) any later version.
 *
 * tiberCAD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with tiberCAD. If not, see <http://www.gnu.org/licenses/>.
 */

/*!
 * \file LatticeMismatch.C
 * \brief tiberCAD elasticity module implementation.
 *
 * \note This file is part of module elasticity.
 */


#include "LatticeMismatch.h"
#include "tibercad/atomistic/BulkCrystal.h"
#include "tibercad/physics/Material.h"
#include "tibercad/module/SimulationInterface.h"
#include "tibercad/module/SimulationEnvironment.h"
#include "tibercad/io/Messages.h"
#include "tibercad/math/TensorOperators.h"

#include "tibercad/module/TiberModule.h"

#include "libmesh/tensor_value.h"


using namespace std;


LatticeMismatch::LatticeMismatch(const ModelOptions& options) :
    BodyForceModel(options)
{
}

LatticeMismatch::~LatticeMismatch(void)
{
  delete _common_supercell;
}


void
LatticeMismatch::do_init(void)
{


  libMesh::RealGradient body_force(0);
  
  //Get reference lattice
  std::string name;
  ModelOptions matopts;

  if (get_options().has_submodel("reference_material"))
  {
    auto ref = get_options().submodels_begin("reference_material");
    matopts = ref->second;

    if (!matopts.find_option("material"))
    {
      matopts["material"] = matopts.get_name();
    }
    
    name = matopts["material"];
  }
  else
  {
    name = get_option("reference_material", "");

    SimulationInterface *sim = SimulationInterface::get_simulation(get_simulator_id());
    const Device& dev = sim->get_environment().get_device();
    const Material* mat = dev.get_material(name);

    if (mat != nullptr)
    {
      matopts = mat->get_options();
      name = mat->get_name();
    }
    else if (!name.empty())
    {
      matopts = get_options();

      get_option("x-growth-direction", "");
      get_option("y-growth-direction", "");
      get_option("z-growth-direction", "");
      get_option("euler_angles", "");
      get_option("x", "");
    }
  }

  if (name.empty())
    throw InitFailedException("Lattice mismatch: reference material is not defined");


  // for "official" materials (created from Device) this is always present
  unsigned int dim = get_material()->get_options().get_option("dimension", 3);
  matopts.set_option("dimension", dim);

  Material* ref_mat = Material::create(name, matopts);
  ref_mat->init();
  const BulkCrystal* ref_bulk = ref_mat->get_bulk_crystal();

  const Material* mat = get_material();
  const BulkCrystal* mat_bulk = mat->get_bulk_crystal();

  libMesh::RealTensor eps0(0);

  if ((mat_bulk != nullptr) && (ref_bulk != nullptr))
  {
    bool compatible = (mat_bulk->get_lattice_type() == ref_bulk->get_lattice_type());

    if (compatible)
    {
      // this one is strained to bulk crystal, and strain tensor rotated
      // to calculation system.

      // cartesian basis vectors in terms of conventional cell are
      // (x y z) = inv(_conv_vec)*_conv_vec
      // so we calculate strained (x' y' z') by using
      // A = inv(_conv_vec) of the strained layer and
      // for _conv_vec the one of the substrate, i.e. we strain the
      // conventinal vectors and calculate the according strain in
      // cartesian coordinates:
      // (x', y', z')' = (1 + eps)*(x, y, z)'
      // => eps = (x', y', z') - 1, because (x, y, z) = 1

      Tensor2 A(inv(mat_bulk->get_conv_vec()));

      Tensor2 eps = A * ref_bulk->get_conv_vec() - Tensor2(1.0);

      // symmetrize to eliminate possible rotations
      eps *= 0.5;
      eps += eps.transpose();

      // rotate to calculation system
      eps = mat_bulk->get_rotation() * eps * mat_bulk->get_rotation().transpose();

      eps0(0, 0) = eps(1, 1);
      eps0(1, 1) = eps(2, 2);
      eps0(2, 2) = eps(3, 3);
      eps0(1, 2) = eps0(2, 1) = eps(3, 2);
      eps0(0, 2) = eps0(2, 0) = eps(3, 1);
      eps0(0, 1) = eps0(1, 0) = eps(2, 1);
    }
    else
    {
      double residual_strain = get_option("max_residual_strain", 0.0);
      delete _common_supercell;
      _common_supercell = new LatticeMatch::Candidate(
          mat_bulk->find_least_common_lattice(*ref_bulk,
                                              eps0, residual_strain));
    }
  }

  double relax = get_option("relaxation_factor", 1.0);

  eps0 *= relax;


  set_strain_source(eps0);

  delete ref_mat;
 
}


void
LatticeMismatch::do_print_info()
{
  Messages msg;

  std::ostringstream os;
  os << "Lattice mismatch body force model" << std::endl;
  os << "Reference material : " << get_option("reference_material", "") << std::endl;
  Messages::info(os.str());
  os.str("");

  if (_common_supercell != nullptr)
  {
    msg.info("Common supercell : ");
    msg.indent();
    msg.info("substrate:");
    msg.indent();
    os << _common_supercell->M1 << std::endl;
    msg.info(os.str());
    msg.unindent();
    os.str("");
    msg.info("superstrate:");
    msg.indent();
    msg.info(os.str());
    os << _common_supercell->M2 << std::endl;
    msg.info(os.str());
    msg.unindent();
    os.str("");
  }

  if (std::abs(get_strain_source().tr()) > 1e-12)
  {
    double relax = get_option("relaxation_factor", 1.0);
    if (relax != 1.0)
      os << "Relaxation factor : " << relax << std::endl;
    msg.info(os.str());
    os.str("");

    msg.info("lattice matching strain:");
    msg.indent();
    os << get_strain_source() << std::endl;
    msg.info(os.str());
    msg.unindent();
  }

}


 
