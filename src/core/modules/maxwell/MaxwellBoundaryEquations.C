// $Id: MaxwellBoundaryEquations.C 2063 2010-09-03 13:11:49Z maufder $

#include "MaxwellBoundaryEquations.h"
#include "SimulationEnvironment.h"
#include "Material.h"
#include "TiberMath.h"
#include "Database.h"
#include "MaxwellBoundaryProperties.h"
#include "OpticPropsModel.h"

#include "equation_systems.h"
#include "dense_submatrix.h"
#include "limits.h"
#include "Utils.h"
#include "TiberLinearSystem.h"
#include "VectorFEBase1D.h"
#include "VectorLinearSystem.h"
#include "OpticPropsInterface.h"

#include "TiberModule.h"

using namespace std;
using namespace Constants;

MaxwellBoundaryEquations::MaxwellBoundaryEquations(const ModelOptions& options)
 : MaxwellEquationsCommon(options)
{
  has_solution_vector(false);
}

//=====================================================//
BoundaryProperties* MaxwellBoundaryEquations::create_boundary_model(const ModelOptions& options) const  throw (ModelErrorException)
{
  return MaxwellBoundaryProperties::create(options);
}

//=======================================================================================================//

void MaxwellBoundaryEquations::do_init() {
  approxOrder = get_options().get_option("approxOrder", 0);
  extraQOrder = get_options().get_option("extraQOrder", 1);
  inplane = get_options().get_option("inplane", "yes");

  W = get_options().get_option("W", 0.0) / Constants::hbar * Constants::e;

  libMesh::EquationSystems& equation_systems = get_equation_systems();

  equation_systems.add_system<VectorLinearSystem> ("MaxwellBoundary");

  VectorLinearSystem& system = equation_systems.get_system<VectorLinearSystem> ("MaxwellBoundary");
  system.simulationInterface = this;

  MeshBase& mesh = get_mesh();

  defaultSourceDirection = 0;
  if (mesh.mesh_dimension() == 3) {
    system.addVariable(approxOrder, true, extraQOrder, false);
    system.addVariable(approxOrder, false);
  } else if (mesh.mesh_dimension() == 2) {
    system.addVariable(approxOrder, true, extraQOrder, inplane == "yes");

    if (inplane == "yes") {
      system.addVariable(approxOrder, false);
    } else {
      defaultSourceDirection = 2; // z
    }
  } else {
    system.addVariable(approxOrder, true, extraQOrder, false);
    defaultSourceDirection = 1; // y
  }

  system.attach_assemble_function(assemble_maxwell_equations);

  system.init();
}


void MaxwellBoundaryEquations::do_reinit() {
  W = get_options().get_option("W", 0.0) / Constants::hbar * Constants::e;
  //std::cout << "WW = " << W << "\n";
}

//=======================================================================================================//
void MaxwellBoundaryEquations::do_solve() {
  std::ostringstream os;
  //os << "Lambda is " << 2 * M_PI * Constants::c / W;
  Messages::info(os.str());

  libMesh::EquationSystems& equation_systems = get_equation_systems();

  VectorLinearSystem& system = equation_systems.get_system<VectorLinearSystem> ("MaxwellBoundary");

  system.solve();
  system.get_solution(edgeSolution);
}

void
MaxwellBoundaryEquations::do_setup_solution_variables(void)
{
  declare_solution(Efield, VECTOR, NODES, "abs");
  declare_solution(Efield_real, VECTOR, NODES, "abs");
  declare_solution(Efield_imag, VECTOR, NODES, "abs");
  declare_solution(Epsilon,
      SolutionDescriptor::TENSOR, SolutionDescriptor::NODES, "abs");
  declare_solution(Epsilon_imag,
      SolutionDescriptor::TENSOR, SolutionDescriptor::NODES, "abs");
  declare_solution(Mu,
      SolutionDescriptor::REAL, SolutionDescriptor::NODES, "abs");
  declare_solution(SVector,
      SolutionDescriptor::VECTOR, SolutionDescriptor::NODES, "abs");

  declare_solution(Bfield,
      SolutionDescriptor::VECTOR, SolutionDescriptor::NODES, "abs");
  declare_solution(Bfield_real,
      SolutionDescriptor::VECTOR, SolutionDescriptor::NODES, "abs");
  declare_solution(Bfield_imag,
      SolutionDescriptor::VECTOR, SolutionDescriptor::NODES, "abs");
  declare_solution(Poynting,
      SolutionDescriptor::VECTOR, SolutionDescriptor::NODES, "abs");
}

// TODO: remove code duplication
void
MaxwellBoundaryEquations::assemble_maxwell_equations(libMesh::EquationSystems& es,
    const std::string& system_name)
{

  const MeshBase& mesh = es.get_mesh();

  const unsigned int dimension = mesh.mesh_dimension();

  VectorLinearSystem& system = es.get_system<VectorLinearSystem>("MaxwellBoundary");

  MaxwellBoundaryEquations* simulation = dynamic_cast<MaxwellBoundaryEquations*>(system.simulationInterface);

  double W = simulation->W;
  double total_length_scaling = system.simulationInterface->get_environment().get_device().get_mesh_units() * system.getGeometryEx()->getScaling().get_length_scaling();
  double K = W / Constants::c * total_length_scaling;

  EdgeDofMap* dof_map = system.getEdgeDofMap(false);

  const VariableType fe_type = system.getVariableType(0);

  IVectorFEBase* fe = dynamic_cast<IVectorFEBase*>(system.getVariableType(0).getFEbase());
  IVectorFEBase* fe_face = dynamic_cast<IVectorFEBase*>(system.getVariableType(0).getFEbase());

  //TODO to be revised
  libMesh::UniquePtr<libMesh::QBase> qrule(new libMesh::QGauss(dimension, static_cast<Order>(2 * fe_type.order + 2 + fe_type.extraQOrder)));

  libMesh::UniquePtr<libMesh::QBase> qrule_face(new libMesh::QGauss(dimension - 1, static_cast<Order>(2 * fe_type.order + 2 + fe_type.extraQOrder)));

  IScalarFEBase* fe_sc = NULL;
  if (system.getVariablesCount() > 1) {
    fe_sc = dynamic_cast<IScalarFEBase*>(system.getVariableType(1).getFEbase());
  }

  MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

  PML pml = system.getGeometryEx()->pml;

  //addBValue(Complex(1, 0), (systemSize - 1) / 2);


  for (; el != end_el; ++el) {
    //std::cout << "Elem\n";
    const Elem* elem = *el;

    std::vector<unsigned int> all_dof_indices;
    system.dof_indices (elem, all_dof_indices);

    fe->attach_quadrature_rule (qrule.get());
    fe_face->attach_quadrature_rule (qrule_face.get());

    if (fe_sc != NULL) {
      fe_sc->attach_quadrature_rule      (qrule.get());
    }

    fe->reinit (elem, dof_map->getPOrder(elem, 0));
    if (fe_sc != NULL) {
      fe_sc->reinit(elem, dof_map->getPOrder(elem, 1));
    }
    const std::vector<Real>& JxW = fe->get_JxW();
    const std::vector<VectorFunction >& edge_phi = fe->getFunctions();
    const std::vector<Point>& xyz = fe->get_xyz();
    const std::vector<Point>& xyz_face = fe_face->get_xyz();

    OpticPropsInterface* opticModel = simulation->getOpticModel(elem);

    //This part is the slowest in assembling.
    for (unsigned int i=0; i<edge_phi.size(); i++) {
      if (all_dof_indices[i] != ElementUtils::INVALID_FUNCTION_ID) {
        for (unsigned int j=0; j<edge_phi.size(); j++) {
          if (all_dof_indices[j] != ElementUtils::INVALID_FUNCTION_ID) {
            libMesh::Complex aValue = 0;

            for (unsigned int qp=0; qp<qrule->n_points(); qp++) {
              libMesh::Complex sInvertDet = pml.getSVectorDet(xyz[qp], opticModel->get_spml());
              aValue += 1/opticModel->get_permeability_constant() * JxW[qp] * (pml.curls(edge_phi[i], xyz[qp], qp, opticModel->get_spml()) * pml.curls(edge_phi[j], xyz[qp], qp, opticModel->get_spml())) / sInvertDet;

              aValue -= simulation->multiply(edge_phi[i].phi[qp], edge_phi[j].phi[qp], opticModel->get_optical_epsilon()) * JxW[qp] / sInvertDet * K * K;
              //aValue -= JxW[qp]*(edge_phi[i].phi[qp] * edge_phi[j].phi[qp])*opticModel->get_dielectric_constant() / sInvertDet * K * K;
            }

            system.addAValue(aValue, all_dof_indices[i], all_dof_indices[j]);
          }
        }
      }

      if (fe_sc != NULL) {
        const std::vector<ScalarFunction >& scalar_phi = fe_sc->getFunctions();

        for (unsigned int j = 0; j < scalar_phi.size(); j++) {
          if (all_dof_indices[i] != ElementUtils::INVALID_FUNCTION_ID && all_dof_indices[j + edge_phi.size()] != ElementUtils::INVALID_FUNCTION_ID) {

            libMesh::Complex aValue = 0;
            for (unsigned int qp=0; qp<qrule->n_points(); qp++) {
              aValue += simulation->multiply(edge_phi[i].phi[qp], scalar_phi[j].grads(qp, pml.getSVector(xyz[qp], opticModel->get_spml())),  opticModel->get_optical_epsilon()) * JxW[qp];
            }

            system.addAValue(aValue, all_dof_indices[i], all_dof_indices[j + edge_phi.size()]);
            system.addAValue(aValue, all_dof_indices[j + edge_phi.size()], all_dof_indices[i]);
          }
        }
      }

      //Apply BC
      for (unsigned int side = 0; side < elem->n_sides(); side++) {
        const libMesh::UniquePtr<libMesh::Elem> sideElem = elem->build_side(i);
        Boundary* bd = simulation->get_environment().get_boundary(ElementSide(elem,i));

        if (bd != NULL && (bd->get_boundary_properties( simulation->get_id() ) != NULL )) {
          MaxwellBoundaryProperties* boundaryProps = dynamic_cast<MaxwellBoundaryProperties*>(bd->get_boundary_properties(simulation->get_id()));
          if (boundaryProps->isSource()) {
            // Add smth...
            fe_face->reinit(elem, dof_map->getPOrder(elem, 0), side);
            const std::vector<Real>& JxW_face = fe_face->get_JxW();
            const std::vector<VectorFunction >& edge_phi_face = fe_face->getFunctions();

            int direction = (boundaryProps->direction < 0 || boundaryProps->direction > 2) ? simulation->defaultSourceDirection : boundaryProps->direction;

            for (unsigned int i=0; i<edge_phi_face.size(); i++) {
              if (all_dof_indices[i] != ElementUtils::INVALID_FUNCTION_ID) {
                libMesh::Complex value = 0;
                    for (unsigned int qp=0; qp<qrule_face->n_points(); qp++) {
                      libMesh::Complex sInvertDet = pml.getSVectorDet(xyz_face[qp], opticModel->get_spml()), one(1, 0);
                      libMesh::Complex sDet = one / sInvertDet;

                      libMesh::Point sourceVector;
                      sourceVector(direction) = boundaryProps->power * K / 2 * system.getGeometryEx()->getScaling().get_length_scaling();//TODO DIRECTION
                      value += sDet * JxW_face[qp] * (edge_phi_face[i].phi[qp] * sourceVector);
                    }
                    system.addBValue(value, all_dof_indices[i]);
              }
            }
          }
        }
      }
    }
  }
  delete fe;
  delete fe_face;
  if (fe_sc != NULL) {
    delete fe_sc;
  }
}

//TODo code duplication
void MaxwellBoundaryEquations::get_solution_secure(const Elem* elem,
    std::map<ID, std::vector<double> >& solutions,
    const std::vector<libMesh::Point>& points)
{
  //Utils::Timer tt;
  libMesh::EquationSystems& equation_systems = get_equation_systems();
  VectorLinearSystem& system = equation_systems.get_system<VectorLinearSystem> ("MaxwellBoundary");

  std::vector<unsigned int> edge_dof_indices;
  IVectorFEBase* fe = dynamic_cast<IVectorFEBase*>(system.getVariableType(0).getFEbase());
  system.getEdgeDofMap(false)->dof_indices (elem, edge_dof_indices, 0);
  fe->reinit(elem, system.getEdgeDofMap(false)->getPOrder(elem, 0), &points);
  const std::vector<VectorFunction>& edge_phi = fe->getFunctions();

  //////////////////////////////

  const std::vector<Point>& xyz = fe->get_xyz();

  PML pml = system.getGeometryEx()->pml;


    for (unsigned int qp = 0; qp < points.size(); qp++) {
      libMesh::Point realValue;
      libMesh::Point imagValue;
      for (unsigned int j = 0; j < edge_dof_indices.size(); j++) {
        if (edge_dof_indices[j] != ElementUtils::INVALID_FUNCTION_ID) {
          realValue +=edgeSolution[edge_dof_indices[j]].real() * edge_phi[j].phi[qp];
          imagValue +=edgeSolution[edge_dof_indices[j]].imag() * edge_phi[j].phi[qp];
        }
      }
  OpticPropsInterface* opticModel = getOpticModel(elem);
/*
    for (unsigned int qp = 0; qp < points.size(); qp++) {
      libMesh::Point realValue;
      libMesh::Point imagValue;
      for (unsigned int j = 0; j < edge_dof_indices.size(); j++) {
        if (edge_dof_indices[j] != ElementUtils::INVALID_FUNCTION_ID) {
          realValue +=edgeSolution[edge_dof_indices[j]].real() * edge_phi[j].phi[qp];
          imagValue +=edgeSolution[edge_dof_indices[j]].imag() * edge_phi[j].phi[qp];
        }
      }
*/

  for (unsigned int qp = 0; qp < points.size(); qp++) {
    if (solutions.count(Mu)) {
      solutions[Mu][qp] = opticModel->get_permeability_constant();
    }

    addTensorSolutionR(solutions, Epsilon, opticModel->get_optical_epsilon(), 6*qp);
    addTensorSolutionI(solutions, Epsilon_imag, opticModel->get_optical_epsilon(), 6*qp);

    if (solutions.count(SVector)) {
      Complex one(1, 0);
      VectorValue<Complex> sVector = pml.getSVector(xyz[qp], opticModel->get_spml());

      solutions[SVector][qp*3] = (one / sVector(0)).imag();
      solutions[SVector][qp*3 + 1] = (one / sVector(1)).imag();
      solutions[SVector][qp*3 + 2] = (one / sVector(2)).imag();
    }
  }

  /////////////////////////////


  for (unsigned int qp = 0; qp < points.size(); qp++) {
    VectorValue<Complex> E_value;
    VectorValue<Complex> B_value;

    for (unsigned int j = 0; j < edge_dof_indices.size(); j++) {
      if (edge_dof_indices[j] != ElementUtils::INVALID_FUNCTION_ID) {
        E_value += getVectorValue(edge_phi[j].phi[qp]) * edgeSolution[edge_dof_indices[j]];
        B_value += edge_phi[j].curl(qp) * edgeSolution[edge_dof_indices[j]];
      }
    }

    B_value = B_value / Complex(0, 1) / W / system.simulationInterface->get_environment().get_device().get_mesh_units() / system.getGeometryEx()->getScaling().get_length_scaling();

      if (solutions.count(SVector)) {
        libMesh::Complex one(1, 0);
        libMesh::VectorValue<libMesh::Complex> sVector = pml.getSVector(xyz[qp], opticModel->get_spml());

    addVectorSolutionR(solutions, Efield_real, E_value, qp*3);
    addVectorSolutionI(solutions, Efield_imag, E_value, qp*3);
    addVectorSolutionA(solutions, Efield, E_value, qp*3);

    addVectorSolutionR(solutions, Bfield_real, B_value, qp*3);
    addVectorSolutionI(solutions, Bfield_imag, B_value, qp*3);
    addVectorSolutionA(solutions, Bfield, B_value, qp*3);

    double mu0 = 4 * M_PI * 1e-7; //Si units
    VectorValue<Complex> H_value = B_value * (1 / opticModel->get_permeability_constant() / mu0);

    VectorValue<Complex> H_conj(std::conj(H_value(0)), std::conj(H_value(1)), std::conj(H_value(2)));
    VectorValue<Complex> P_value = E_value.cross(H_conj) * 0.5;

    addVectorSolutionR(solutions, Poynting, P_value, qp*3);
  }

  delete fe;
}

OpticPropsInterface* MaxwellBoundaryEquations::getOpticModel(const Elem* elem) {
  ID subdomain = elem->subdomain_id();
  const Material* material = get_environment().get_device().get_material(subdomain);

  OpticPropsInterface* result = dynamic_cast<OpticPropsInterface*>(
          material->get_model(get_id()));


  return result;
}

