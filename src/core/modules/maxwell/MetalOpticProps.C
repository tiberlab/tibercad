
#include "MetalOpticProps.h"
#include "SimulationInterface.h"
#include "Material.h"
#include "TiberModule.h"

void MetalOpticProps::read_database() {
  OpticPropsModel::read_database();

  const Database& db = get_database();
  db.set_section("permittivity");

  double alpha_real = db.get("optical_alpha", 0.0);

  double alpha_imag = db.get("optical_alpha_imag", 0.0);

  double W1_real = db.get("optical_W1", 0.0);

  double W1_imag = db.get("optical_W1_imag", 0.0);

  alpha = Complex(alpha_real, alpha_imag);

  W1 = Complex(W1_real, W1_imag);
}

void MetalOpticProps::do_init(void) {
  get_parameter_c("optical_alpha", alpha);
  get_parameter_c("optical_W1", W1);
  OpticPropsModel::do_init();
  MaxwellEquations* maxwell = dynamic_cast<MaxwellEquations*>(SimulationInterface::get_simulation(get_simulator_id()));
  maxwell->cubics.push_back(*this);
}

bool MetalOpticProps::addCData(CubicEigenSystem& system) {


  MaxwellEquations* maxwell = dynamic_cast<MaxwellEquations*>(SimulationInterface::get_simulation(get_simulator_id()));

  const MeshBase& mesh = system.get_mesh();
  const unsigned int dimension = mesh.mesh_dimension();

  IVectorFEBase* fe = dynamic_cast<IVectorFEBase*>(system.getVariableType(0).getFEbase());
  const VariableType fe_type = system.getVariableType(0);
  //TODO to be revised
  UniquePtr<QBase> qrule(new QGauss(dimension, static_cast<Order>(2 * fe_type.order + 2 + fe_type.extraQOrder)));

  MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

  for (; el != end_el; ++el) {
    const Elem* elem = *el;

    PML pml = system.getGeometryEx()->pml;

    OpticPropsInterface* opticModel = maxwell->getOpticModel(elem);

    if (!pml.isPMLRegion(elem, maxwell)) {

      std::vector<unsigned int> all_dof_indices;
      system.dof_indices (elem, all_dof_indices);

      fe->attach_quadrature_rule (qrule.get());
      fe->reinit (elem, system.getEdgeDofMap(false)->getPOrder(elem, 0));

      const std::vector<Real>& JxW = fe->get_JxW();
      const std::vector<VectorFunction >& edge_phi = fe->getFunctions();
      const std::vector<Point>& xyz = fe->get_xyz();

      for (unsigned int i=0; i<edge_phi.size(); i++) {
        if (all_dof_indices[i] != ElementUtils::INVALID_FUNCTION_ID) {
          for (unsigned int j=0; j<edge_phi.size(); j++) {
            if (all_dof_indices[j] != ElementUtils::INVALID_FUNCTION_ID) {
              Complex bValue = 0;

              for (unsigned int qp=0; qp<qrule->n_points(); qp++) {
                Complex sInvertDet = pml.getSVectorDet(xyz[qp], opticModel->get_spml());
                bValue += maxwell->multiply(edge_phi[i].phi[qp], edge_phi[j].phi[qp], opticModel->get_optical_epsilon()) * JxW[qp] / sInvertDet;
              }

              system.addCValue(bValue, all_dof_indices[i], all_dof_indices[j]);
            }
          }
        }
      }
    }
  }

  delete fe;
}
