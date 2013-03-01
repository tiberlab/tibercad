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
 : SimulationInterface(options)
{
  has_solution_vector(false);
}

//=====================================================//
BoundaryProperties* MaxwellBoundaryEquations::create_boundary_model(const ModelOptions& options) const  throw (ModelErrorException)
{
  return MaxwellBoundaryProperties::create(options);
}

//=======================================================================================================//
PhysicalModel*  MaxwellBoundaryEquations::create_physical_model(const ModelOptions& options,
    const Material* mat) const throw (ModelErrorException)
{
  const std::string& modelname = get_option((mat->get_name() + "_opticmodel"), "");

  if (modelname == "") {
    return OpticPropsModel::create(options);
  } else {
    OpticPropsInterface* model =
        OpticPropsInterface::create(modelname, mat, options);

    if (model == NULL)
      throw ModelErrorException(
          "Maxwell: No such physical model: " + modelname);

    return model;
  }
}

//=======================================================================================================//

void MaxwellBoundaryEquations::do_init() {
  approxOrder = get_options().get_option("approxOrder", 0);
  extraQOrder = get_options().get_option("extraQOrder", 1);
  inplane = get_options().get_option("inplane", "yes");

  W = get_options().get_option("W", 0.0) / Constants::hbar * Constants::e;

  EquationSystems& equation_systems = get_equation_systems();

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
}

//=======================================================================================================//
void MaxwellBoundaryEquations::do_solve() {
  std::ostringstream os;
  os << "Lambda is " << 2 * M_PI * Constants::c / W;
  Messages::info(os.str());

  EquationSystems& equation_systems = get_equation_systems();

  VectorLinearSystem& system = equation_systems.get_system<VectorLinearSystem> ("MaxwellBoundary");

  system.solve();
  system.get_solution(edgeSolution);
}

void
MaxwellBoundaryEquations::do_setup_solution_variables(void)
{
  declare_solution(Intensity, SolutionDescriptor::REAL, NODES, "abs");
  declare_solution(Efield, VECTOR, NODES, "abs");
  declare_solution(Efield_real, VECTOR, NODES, "abs");
  declare_solution(Efield_imag, VECTOR, NODES, "abs");
  declare_solution(Epsilon,
      SolutionDescriptor::REAL, SolutionDescriptor::NODES, "abs");
  declare_solution(Epsilon_imag,
      SolutionDescriptor::REAL, SolutionDescriptor::NODES, "abs");
  declare_solution(Mu,
      SolutionDescriptor::REAL, SolutionDescriptor::NODES, "abs");
  declare_solution(SVector,
      SolutionDescriptor::VECTOR, SolutionDescriptor::NODES, "abs");
}

// TODO: remove code duplication
void
MaxwellBoundaryEquations::assemble_maxwell_equations(EquationSystems& es,
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
  AutoPtr<QBase> qrule(new QGauss(dimension, static_cast<Order>(2 * fe_type.order + 2 + fe_type.extraQOrder)));

  AutoPtr<QBase> qrule_face(new QGauss(dimension - 1, static_cast<Order>(2 * fe_type.order + 2 + fe_type.extraQOrder)));

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
            Complex aValue = 0;
            //Complex bValue = 0;

            for (unsigned int qp=0; qp<qrule->n_points(); qp++) {
              Complex sInvertDet = pml.getSVectorDet(xyz[qp], opticModel->get_spml());
              //std::cout << "spml: " << params.sPML << "\n";
              aValue += 1/opticModel->get_permeability_constant() * JxW[qp] * (pml.curls(edge_phi[i], xyz[qp], qp, opticModel->get_spml()) * pml.curls(edge_phi[j], xyz[qp], qp, opticModel->get_spml())) / sInvertDet;

              aValue -= simulation->multiply(edge_phi[i].phi[qp], edge_phi[j].phi[qp], opticModel->get_optical_epsilon(), opticModel->get_optical_epsilon_imag()) * JxW[qp] / sInvertDet * K * K;
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

            Complex aValue = 0;
            for (unsigned int qp=0; qp<qrule->n_points(); qp++) {
              aValue += simulation->multiply(edge_phi[i].phi[qp], scalar_phi[j].grads(qp, pml.getSVector(xyz[qp], opticModel->get_spml())(0)),  opticModel->get_optical_epsilon(), opticModel->get_optical_epsilon_imag()) * JxW[qp];
            }

            system.addAValue(aValue, all_dof_indices[i], all_dof_indices[j + edge_phi.size()]);
            system.addAValue(aValue, all_dof_indices[j + edge_phi.size()], all_dof_indices[i]);
          }
        }
      }

      //Apply BC
      for (unsigned int side = 0; side < elem->n_sides(); side++) {
        const AutoPtr<Elem> sideElem = elem->build_side(i);
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
                    Complex value = 0;
                    for (unsigned int qp=0; qp<qrule_face->n_points(); qp++) {
                      Complex sInvertDet = pml.getSVectorDet(xyz_face[qp], opticModel->get_spml()), one(1, 0);
                      Complex sDet = one / sInvertDet;

                      Point sourceVector;
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
    const std::vector<Point>& points)
{
  //Utils::Timer tt;
  EquationSystems& equation_systems = get_equation_systems();
  VectorLinearSystem& system = equation_systems.get_system<VectorLinearSystem> ("MaxwellBoundary");

  std::vector<unsigned int> edge_dof_indices;
  IVectorFEBase* fe = dynamic_cast<IVectorFEBase*>(system.getVariableType(0).getFEbase());
  system.getEdgeDofMap(false)->dof_indices (elem, edge_dof_indices, 0);
  fe->reinit(elem, system.getEdgeDofMap(false)->getPOrder(elem, 0), &points);
  const std::vector<VectorFunction>& edge_phi = fe->getFunctions();

 // std::cout << "1:" << tt.elapsed_string() << "\n";
  if (solutions.count(Efield)) {
    std::vector<double>& solution = solutions[Efield];
    std::vector<double>& solution_real = solutions[Efield_real];
    std::vector<double>& solution_imag = solutions[Efield_imag];
    std::vector<double>& solution_intensity = solutions[Intensity];

    //std::vector<Complex> edgeSolution;
    //system.get_solution(edgeSolution);
    //std::cout << "2:" << tt.elapsed_string() << "\n";
/*
    if (edge_dof_indices[0] == 0) {
      //std::cout << "EDGE SOLUTION: \n";
      for (int iii = 0; iii < edgeSolution.size(); iii++) {
        std::cout << iii << " " << edgeSolution[iii];
      }
    }
*/

    solution.resize(points.size()*3);
    solution_real.resize(points.size()*3);
    solution_imag.resize(points.size()*3);

    for (unsigned int qp = 0; qp < points.size(); qp++) {
      Point realValue;
      Point imagValue;
      for (unsigned int j = 0; j < edge_dof_indices.size(); j++) {
        if (edge_dof_indices[j] != ElementUtils::INVALID_FUNCTION_ID) {
          realValue +=edgeSolution[edge_dof_indices[j]].real() * edge_phi[j].phi[qp];
          imagValue +=edgeSolution[edge_dof_indices[j]].imag() * edge_phi[j].phi[qp];
        }
      }

      solution[qp*3] = std::sqrt(realValue(0) * realValue(0) + imagValue(0) * imagValue(0));
      solution[qp*3 + 1] = std::sqrt(realValue(1) * realValue(1) + imagValue(1) * imagValue(1));
      solution[qp*3 + 2] = std::sqrt(realValue(2) * realValue(2) + imagValue(2) * imagValue(2));
//      solution[qp*3] = (realValue(0) * realValue(0) + imagValue(0) * imagValue(0));
//      solution[qp*3 + 1] = (realValue(1) * realValue(1) + imagValue(1) * imagValue(1));
//      solution[qp*3 + 2] = (realValue(2) * realValue(2) + imagValue(2) * imagValue(2));

      solution_real[qp*3] = realValue(0);
      solution_real[qp*3 + 1] = realValue(1);
      solution_real[qp*3 + 2] = realValue(2);

      solution_imag[qp*3] = imagValue(0);
      solution_imag[qp*3 + 1] = imagValue(1);
      solution_imag[qp*3 + 2] = imagValue(2);

      solution_intensity[qp] = realValue * realValue + imagValue * imagValue;
    }
  }
  //std::cout << "3:" << tt.elapsed_string() << "\n";
  if (solutions.count(Mu) || solutions.count(Epsilon) || solutions.count(SVector)) {
    const std::vector<Point>& xyz = fe->get_xyz();

    PML pml = system.getGeometryEx()->pml;

    OpticPropsInterface* opticModel = getOpticModel(elem);

    for (unsigned int qp = 0; qp < points.size(); qp++) {
      if (solutions.count(Mu)) {
        solutions[Mu][qp] = opticModel->get_permeability_constant();
      }

      if (solutions.count(Epsilon)) {
        solutions[Epsilon][qp] = opticModel->get_dielectric_constant().real();
      }

      if (solutions.count(Epsilon_imag)) {
        solutions[Epsilon_imag][qp] = opticModel->get_dielectric_constant().imag();
      }

      if (solutions.count(SVector)) {
        Complex one(1, 0);
        VectorValue<Complex> sVector = pml.getSVector(xyz[qp], opticModel->get_spml());

        solutions[SVector][qp*3] = (one / sVector(0)).imag();
        solutions[SVector][qp*3 + 1] = (one / sVector(1)).imag();
        solutions[SVector][qp*3 + 2] = (one / sVector(2)).imag();
      }
    }
  }

  delete fe;
}

OpticPropsInterface* MaxwellBoundaryEquations::getOpticModel(const Elem* elem) {
  ID subdomain = elem->subdomain_id();
  const Material* material = get_environment().get_device().get_material(subdomain);

  OpticPropsInterface* result = dynamic_cast<OpticPropsInterface*>(
          material->get_model(get_id()));

  return result;

  return result;
}

Complex MaxwellBoundaryEquations::multiply(const Point& v1, const Point& v2, const Tensor2Sym& t_real, const Tensor2Sym& t_imag) {
  Complex result = 0;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      Complex t(t_real(i+1, j+1), t_imag(i+1, j+1));

      result += t * v1(i) * v2(j);
    }
  }

  return result;
}

Complex MaxwellBoundaryEquations::multiply(const Point& v1, const VectorValue<Complex>& v2, const Tensor2Sym& t_real, const Tensor2Sym& t_imag) {
  Complex result = 0;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      Complex t(t_real(i+1, j+1), t_imag(i+1, j+1));
      result += t * v1(i) * v2(j);
    }
  }

  return result;
}
