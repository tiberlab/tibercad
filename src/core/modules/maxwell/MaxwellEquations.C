// $Id: MaxwellEquations.C 2063 2010-09-03 13:11:49Z maufder $

#include "MaxwellEquations.h"
#include "EigenSolver.h"
#include "SimulationEnvironment.h"
#include "Material.h"
#include "TiberMath.h"
#include "Database.h"
#include <assert.h>
#include "MaxwellBoundaryProperties.h"
#include "OpticPropsModel.h"

#include "equation_systems.h"
#include "dense_submatrix.h"
#include "VectorFunction.h"
#include "ScalarFunction.h"
#include "MaxwellBoundaryProperties.h"
#include "IVectorFEBase.h"
#include "IScalarFEBase.h"
#include "limits.h"
#include "Utils.h"
using namespace std;
using namespace Constants;
#include "sys/time.h"
#include "ExcitonLayer.h"
#include "PMLFilter.h"


//TODO
//TODO add tensor here
//TODO remove code duplication
//TODO

TIBER_MODULE(MaxwellEquations, MODULE_NAME)

MaxwellEquations::MaxwellEquations(const ModelOptions& options)
 : SimulationInterface(options)
{
  has_solution_vector(false);
  Wc0 = Complex(-1, 0);
}

//=====================================================//
BoundaryProperties* MaxwellEquations::create_boundary_model(const ModelOptions& options) const  throw (ModelErrorException)
{
  return MaxwellBoundaryProperties::create(options);
}

//=======================================================================================================//
PhysicalModel*  MaxwellEquations::create_physical_model(const ModelOptions& options,
    const Material* mat) const throw (ModelErrorException)
{
  return OpticPropsModel::create(options);
}

//=======================================================================================================//

void MaxwellEquations::do_init() {
  eigenCount = get_options().get_option("eigenCount", 15);

  approxOrder = get_options().get_option("approxOrder", 0);
  extraQOrder = get_options().get_option("extraQOrder", 1);

  maxIterations = get_options().get_option("maxIterations", 1);
  relativeError = get_options().get_option("relativeError", 0.2);
  inplane = get_options().get_option("inplane", "yes");

  spectrumShift = get_options().get_option("spectrumShift", 0.0);

  solver_max_it = get_options().get_option("solver_maxits", 100);
  solver_tolerance = get_options().get_option("solver_rtol", 0.001);
  storeSolutions = get_options().get_option("storeSolutions", eigenCount < 100);

  polaritons = get_options().get_option("polaritons", false);
  wellMaterial = get_options().get_option("wellMaterial", "GaN");
  lwell = get_options().get_option("lwell", 3.0) * 1e-9;
  errorEstamite = get_options().get_option("errorEstimate", false);
  useCubic = false;

  Wlt = get_options().get_option("Wlt", 0.0) / Constants::hbar * Constants::e;

  double Wexc0_r = get_options().get_option("Wexc", 0.0) / Constants::hbar * Constants::e;
  double Wexc0_i = - 1.0 / get_options().get_option("excitonLifeTime", 1.0);
  Wexc0 = Complex(Wexc0_r, Wexc0_i);

  MeshBase& mesh = get_mesh();

  EquationSystems& equation_systems = get_equation_systems();


  if (polaritons && useCubic) {
    equation_systems.add_system<CubicEigenSystem> ("Maxwell");
  } else {
    equation_systems.add_system<EigenSystem> ("Maxwell");
  }

  EigenSystem& system = equation_systems.get_system<EigenSystem> ("Maxwell");
  system.simulationInterface = this;

  system.setSpectrumShift(Complex(spectrumShift * system.simulationInterface->get_environment().get_device().get_mesh_units() / c, 0));

  system.solver_max_it = solver_max_it;
  system.solver_tolerance = solver_tolerance;

  //std::cout << "Spectrum shift: " << system.spectrumShift << "\n";

  system.setRequestedEigenPairs(eigenCount);

  if (mesh.mesh_dimension() == 3) {
    system.addVariable(approxOrder, true, extraQOrder, false);
    system.addVariable(approxOrder, false);
  } else if (mesh.mesh_dimension() == 2) {
    system.addVariable(approxOrder, true, extraQOrder, inplane == "yes");

    if (inplane == "yes") {
      system.addVariable(approxOrder, false);
    }
  } else {
    system.addVariable(approxOrder, true, extraQOrder, false);
  }

  system.attach_assemble_function(assemble_maxwell_equations);
  relativeIndexing = false;
}


//=======================================================================================================//
void MaxwellEquations::do_solve() {
  MeshBase& mesh = get_mesh();

  EquationSystems& equation_systems = get_equation_systems();

  EigenSystem& system = equation_systems.get_system<EigenSystem> ("Maxwell");

  if (polaritons && !useCubic) { // We already obtained Wc0 & Some hopfield coefficients
    if (Wc0.real() > 0) {
/*
      ExcitonLayer layer(system.getEdgeDofMap(false)->getGeometryEx()->getScaling().get_length_scaling() *
          system.simulationInterface->get_environment().get_device().get_mesh_units(), this);

      Complex newVRabiApprox = layer.getApproximateRabiSplitting(system);
      if (VRabiApprox.real() > 0) {
        std::cout << "NEW RABI APPROX IN eV: " <<  newVRabiApprox * Constants::hbar / Constants::e << "\n";
        std::cout << "OLD RABI APPROX IN eV: " <<  VRabiApprox * Constants::hbar / Constants::e << "\n";
        std::cout << "DIFF IN %: " <<  std::abs((VRabiApprox - newVRabiApprox)/VRabiApprox) << "\n";
        if (std::abs((VRabiApprox - newVRabiApprox)/VRabiApprox) <= 0.01) {
          std::cout << "Polariton mode didnt changed\n";
          return;
        } else {
          VRabiApprox = newVRabiApprox;
          std::cout << "Polariton mode changed\n";
          //return;
        }
        flush(std::cout);
      }
      VRabiApprox = newVRabiApprox;
*/
      calculateHopfieldCoefficients();
      return;
    }
  }



  int stime = time(NULL);

  system.reinit();

  system.solve();

  if (storeSolutions) {
    for (int i = 0; i < system.get_n_converged(); i++) {
      system.get_eigen_vector(i, storedSolutions[i]);
    }
  }
  increment_solve_sequence_number();

  filterEigenValues(0.0);

  int ttime = time(NULL) - stime;
  //cout << "Solving time in seconds: " << ttime << "\n";
  //flush(cout);

  std::ostringstream os;

  std::vector<Complex> eigens;
  for (int i = 0; i < accepted_eigen_count; i++) {
    Complex eigen1 = system.get_eigen_lambda(eigenIndices[i]);
    eigens.push_back(eigen1);

    Complex eigen2 = eigen1 * c / system.simulationInterface->get_environment().get_device().get_mesh_units();

    os << i << " " << eigen1 << " W= || " <<  eigen2 << " || or in eV: " << (eigen2  * Constants::hbar / Constants::eV) << "\n";
  }

  os << "Eigen total: " << system.get_n_converged() << ", eigen accepted: " << accepted_eigen_count << ".\n";
//  cout << "======================================================================\n";
  Messages::info(os.str());

  if (Wc0.real() == -1) {
    Wc0 =  system.get_eigen_lambda(eigenIndices[0])  * c / system.simulationInterface->get_environment().get_device().get_mesh_units();
    if (polaritons) {
      std::cout << "Setting Wc0 to " << Wc0 << "\n";
      flush(std::cout);
    }
  }

  if (errorEstamite) {
    errorEstamite = false;
    approxOrder++;

    EigenSystem& system = equation_systems.get_system<EigenSystem> ("Maxwell");
    for (int i = 0; i < system.variables.size(); i++) {
      system.variables[i].order++;
    }

    do_solve();

    std::ostringstream os1;
    os1 << "Real error can be either bigger or smaller than presented values\n";
    for (int i = 0; i < accepted_eigen_count && i < eigens.size(); i++) {
      Complex eigen1 = system.get_eigen_lambda(eigenIndices[i]);
      double error = std::abs(eigen1 - eigens[i])/std::abs(eigen1);
      os1 << i << " error: " << error << "\n";
    }
    Messages::info(os1.str());
  }
}

void MaxwellEquations::filterEigenValues(double factor) {
  EquationSystems& equation_systems = get_equation_systems();

  EigenSystem& system = equation_systems.get_system<EigenSystem> ("Maxwell");

  PMLFilter filter(factor, this);

  relativeIndexing = false;

  Utils::Timer tt;

  filter.filter(system, eigenIndices, storedSolutions);

  //std::cout << "Filtered in " << tt.elapsed_string() << "\n";
  //flush(std::cout);

  accepted_eigen_count = eigenIndices.size();
  relativeIndexing = true;
}

void
MaxwellEquations::do_setup_solution_variables(void)
{
  eigensOut = get_options().get_option("eigensOut", 1);

  declare_solution(EigenValue, REAL, GLOBAL, "c-1");
  declare_solution(EigenValue_eV, REAL, GLOBAL, "eV");
  declare_solution(EigenValueImag, REAL, GLOBAL, "c-1");

  //if (polaritons) {
    declare_solution(WPolariton, REAL, GLOBAL, "c-1");
    declare_solution(WPolariton_eV, REAL, GLOBAL, "eV");
    declare_solution(WPolaritonImag, REAL, GLOBAL, "c-1");
  //}

  for (int i = 0; i < eigensOut; i++) {
    char buffer [50];
    sprintf(buffer, "Efield_%d", i);

    declare_solution_ext(buffer , Efield + i,
        SolutionDescriptor::VECTOR, SolutionDescriptor::NODES, "AnyLoL");
  }

  declare_solution(XHopfield, SolutionDescriptor::REAL, SolutionDescriptor::NODES, "abs");
  declare_solution(Epsilon, SolutionDescriptor::REAL, SolutionDescriptor::NODES, "abs");
  declare_solution(Epsilon_imag, SolutionDescriptor::REAL, SolutionDescriptor::NODES, "abs");
  declare_solution(Mu, SolutionDescriptor::REAL, SolutionDescriptor::NODES, "abs");
  declare_solution(SVector, SolutionDescriptor::VECTOR, SolutionDescriptor::NODES, "abs");

}



void
MaxwellEquations::assemble_maxwell_equations(EquationSystems& es,
    const std::string& system_name)
{
  const MeshBase& mesh = es.get_mesh();

  const unsigned int dimension = mesh.mesh_dimension();

  EigenSystem& system = es.get_system<EigenSystem>("Maxwell");

  MaxwellEquations* this_mme = dynamic_cast<MaxwellEquations*>(system.simulationInterface);

  // Here we should pass all length scaling
  ExcitonLayer layer(system.getEdgeDofMap(false)->getGeometryEx()->getScaling().get_length_scaling() *
      system.simulationInterface->get_environment().get_device().get_mesh_units(), this_mme);

  if (this_mme->polaritons && this_mme->useCubic) {
    layer.preset(system);
  }

  EdgeDofMap* dof_map = system.getEdgeDofMap(false);

  const VariableType fe_type = system.getVariableType(0);

  IVectorFEBase* fe = dynamic_cast<IVectorFEBase*>(system.getVariableType(0).getFEbase());

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

  for (; el != end_el; ++el) {
    const Elem* elem = *el;

    std::vector<unsigned int> all_dof_indices;
    system.dof_indices (elem, all_dof_indices);

    fe->attach_quadrature_rule (qrule.get());
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

    OpticPropsInterface* opticModel = this_mme->getOpticModel(elem);

    //This part is the slowest in assembling.
    for (unsigned int i=0; i<edge_phi.size(); i++) {
      if (all_dof_indices[i] != ElementUtils::INVALID_FUNCTION_ID) {
        for (unsigned int j=0; j<edge_phi.size(); j++) {
          if (all_dof_indices[j] != ElementUtils::INVALID_FUNCTION_ID) {

            Complex aValue = 0;
            Complex bValue = 0;

            for (unsigned int qp=0; qp<qrule->n_points(); qp++) {
              Complex sInvertDet = pml.getSVectorDet(xyz[qp], opticModel->get_spml());
              aValue += 1/opticModel->get_permeability_constant() * JxW[qp] * (pml.curls(edge_phi[i], xyz[qp], qp, opticModel->get_spml()) * pml.curls(edge_phi[j], xyz[qp], qp, opticModel->get_spml())) / sInvertDet;
              bValue += JxW[qp]*(edge_phi[i].phi[qp] * edge_phi[j].phi[qp])*opticModel->get_dielectric_constant() / sInvertDet;
            }

            system.addAValue(aValue, all_dof_indices[i], all_dof_indices[j]);
            system.addBValue(bValue, all_dof_indices[i], all_dof_indices[j]);
          }
        }
      }

      //TODO add svector here???
      if (fe_sc != NULL) {
        const std::vector<ScalarFunction >& scalar_phi = fe_sc->getFunctions();

        for (unsigned int j = 0; j < scalar_phi.size(); j++) {
          if (all_dof_indices[i] != ElementUtils::INVALID_FUNCTION_ID && all_dof_indices[j + edge_phi.size()] != ElementUtils::INVALID_FUNCTION_ID) {

            Complex aValue = 0;
            for (unsigned int qp=0; qp<qrule->n_points(); qp++) {
              aValue += opticModel->get_dielectric_constant() * JxW[qp] * (edge_phi[i].phi[qp] * scalar_phi[j].grads(qp, pml.getSVector(xyz[qp], opticModel->get_spml())));
            }

            system.addAValue(aValue, all_dof_indices[i], all_dof_indices[j + edge_phi.size()]);
            system.addAValue(aValue, all_dof_indices[j + edge_phi.size()], all_dof_indices[i]);
          }
        }
      }
    }
  }
  delete fe;
  if (fe_sc != NULL) {
    delete fe_sc;
  }


  if (this_mme->polaritons && this_mme->Wc0.real() > 0 && this_mme->useCubic) {
    //std::cout << "Add data Called\n"; flush(std::cout);
    layer.addData(system);
  }
}

void MaxwellEquations::get_solution_secure(const Elem* elem,
    std::map<ID, std::vector<double> >& solutions,
    const std::vector<Point>& points)
{

  EigenSystem& system = get_equation_systems().get_system<EigenSystem>("Maxwell");

  std::vector<unsigned int> edge_dof_indices;
  IVectorFEBase* fe = dynamic_cast<IVectorFEBase*>(system.getVariableType(0).getFEbase());
  system.getEdgeDofMap(false)->dof_indices (elem, edge_dof_indices, 0);
  fe->reinit(elem, system.getEdgeDofMap(false)->getPOrder(elem, 0), &points);
  const std::vector<VectorFunction>& edge_phi = fe->getFunctions();

  int iMax = relativeIndexing ? accepted_eigen_count : system.get_n_converged();

  for (int i = 0; i < iMax; i++) {
    if (solutions.count(Efield + i)) {
      std::vector<double>& solution = solutions[Efield + i];

      int ii = relativeIndexing ? eigenIndices[i] : i;

      std::vector<Complex> edgeSolution;
      if (!storeSolutions) {
        system.get_eigen_vector(ii, edgeSolution);
      }

      solution.resize(points.size()*3);
      for (unsigned int qp = 0; qp < points.size(); qp++) {
        Point value;
        for (unsigned int j = 0; j < edge_dof_indices.size(); j++) {
          if (edge_dof_indices[j] != ElementUtils::INVALID_FUNCTION_ID) {
            if (!storeSolutions) {
              value = value + edgeSolution[edge_dof_indices[j]].real() * edge_phi[j].phi[qp];
            } else {
              value = value + storedSolutions[ii][edge_dof_indices[j]].real() * edge_phi[j].phi[qp];
            }
          }
        }
        solution[qp*3] = value(0);
        solution[qp*3 + 1] = value(1);
        solution[qp*3 + 2] = value(2);
      }
    }
  }

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

  if (solutions.count(XHopfield)) {
    ID subdomain = elem->subdomain_id();
    double hopfield = hopfieldCoeeficients.count(subdomain) ? hopfieldCoeeficients[subdomain] : 0;
    for (unsigned int qp = 0; qp < points.size(); qp++) {
      solutions[XHopfield][qp] = hopfield;
    }
  }

}

void MaxwellEquations::get_solution_secure(std::map<ID, std::vector<double> >& solutions) {
  EigenSystem& system = get_equation_systems().get_system<EigenSystem>("Maxwell");

  if (solutions.count(EigenValue)) {
    std::vector<double>& solution = solutions[EigenValue];
    solution.resize(0);
    for (int i = 0; i < accepted_eigen_count && i < eigensOut; i++) {
      double eigen1 = system.get_eigen_lambda(eigenIndices[i]).real();

      solution.push_back(eigen1 * c / system.simulationInterface->get_environment().get_device().get_mesh_units());
    }
  }

  if (solutions.count(EigenValue_eV)) {
    std::vector<double>& solution = solutions[EigenValue_eV];
    solution.resize(0);
    for (int i = 0; i < accepted_eigen_count && i < eigensOut; i++) {
      double eigen1 = system.get_eigen_lambda(eigenIndices[i]).real();

      solution.push_back(eigen1 * c / system.simulationInterface->get_environment().get_device().get_mesh_units() * Constants::hbar / Constants::e);
    }
  }

  if (solutions.count(EigenValueImag)) {
    std::vector<double>& solution = solutions[EigenValueImag];
    solution.resize(0);
    for (int i = 0; i < accepted_eigen_count && i < eigensOut; i++) {
      double eigen1 = system.get_eigen_lambda(eigenIndices[i]).imag();

      solution.push_back(eigen1 * c / system.simulationInterface->get_environment().get_device().get_mesh_units());
    }
  }

  if (solutions.count(WPolariton)) {
    solutions[WPolariton][0] = WPolaritonLow.real();
  }

  if (solutions.count(WPolariton_eV)) {
    solutions[WPolariton_eV][0] = WPolaritonLow.real() * Constants::hbar / Constants::e;
  }

  if (solutions.count(WPolaritonImag)) {
    solutions[WPolaritonImag][0] = WPolaritonLow.imag();
  }
}

// Return square of module of X0.
/*
double MaxwellEquations::getXHopfield(int i) {
  EquationSystems& equation_systems = get_equation_systems();

  EigenSystem& system = equation_systems.get_system<EigenSystem> ("Maxwell");

  Complex Wlow = system.get_eigen_lambda(eigenIndices[0])  * c / system.simulationInterface->get_environment().get_device().get_mesh_units();


   * We need eigenvectors:
   *
   * Wc  V        C      =  Wlow C
   * V  Wexc0     X      =  Wlow X
   *
   * WcC + VX = WlowC
   * VC + WxX = WlowX
   *
   * WcC + X (WlowX - WxX) / C = WlowC
   * C*C(Wc-Wlow) = X*X*(Wx-Wlow)
   * CC*+XX* = 1
   *
   * XX* (1 + abs(Wx-Wlow)/abs(Wc - Wlow)) = 1

  Complex t(1, 0);
  std::cout << "Hop " << Wc << " " << Wexc0 << " " << Wlow << "\n"; flush(std::cout);

  // V = (WlowC - WcC) / X = (Wlow - Wc) C/X =
  //double V = (Wlow - Wc) * std::sqrt(1 - )
  //std::cout << "OLD RABI " << std::sqrt((Wlow - Wc)*(Wlow - Wexc0)) << "\n"; flush(std::cout);
  std::cout << "PRECISE RABI IN eV: " << std::sqrt((Wlow - Wc)*(Wlow - Wexc0)) * Constants::hbar / Constants::e << "\n";

  return 1.0 / (1.0 + std::abs(((Wexc0 - Wlow)/(Wc - Wlow))));
}
*/

void MaxwellEquations::plot_globaldata() {
  if (!polaritons) {
    SimulationInterface::plot_globaldata();
  }
}

OpticPropsInterface* MaxwellEquations::getOpticModel(const Elem* elem) {
  ID subdomain = elem->subdomain_id();
  const Material* material = get_environment().get_device().get_material(subdomain);


  return dynamic_cast<OpticPropsInterface*>(
          material->get_model(get_id()));

}

void MaxwellEquations::calculateHopfieldCoefficients() {
  SimulationInterface* _exciton_sim = find_simulation("excitontransport");

  if (_exciton_sim == NULL) {
    std::string msg("ExcitonLayer: Simulation not found");
    throw InitFailedException(msg);
  } else if (!_exciton_sim->is_solved() || Wc0.real() == -1) {
    return;
  }

  ID densityId = _exciton_sim->get_solution_id("Xdens");

  double E2 = 0; //square of E abs. (with epsilon weight)
  std::map<ID, double> Fexc2; // square of Fexc abs. For each well.
  std::map<ID, Point> FexcE; // integral of Fexc*E.  For each well.

  EigenSystem& system = get_equation_systems().get_system<EigenSystem>("Maxwell");

  const MeshBase& mesh = system.get_mesh();
  const unsigned int dimension = mesh.mesh_dimension();

  IVectorFEBase* fe = dynamic_cast<IVectorFEBase*>(system.getVariableType(0).getFEbase());
  const VariableType fe_type = system.getVariableType(0);
  //TODO to be revised
  AutoPtr<QBase> qrule(new QGauss(dimension, static_cast<Order>(2 * fe_type.order + 2 + fe_type.extraQOrder)));

  MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

  PML pml = system.getGeometryEx()->pml;

  for (; el != end_el; ++el) {
    const Elem* elem = *el;

    ID subdomain = elem->subdomain_id();
    if (!Fexc2.count(subdomain)) {
      Fexc2[subdomain] = 0;
      FexcE[subdomain] = Point(0);
    }

    //OpticPropsInterface* opticModel = getOpticModel(elem);

    //TODO
    bool inWell = get_environment().get_device().get_material(elem->subdomain_id())->get_name() == wellMaterial;
    if (!pml.isPMLRegion(elem, this)) {

      std::vector<unsigned int> all_dof_indices;
      system.dof_indices (elem, all_dof_indices);

      fe->attach_quadrature_rule (qrule.get());
      fe->reinit (elem, system.getEdgeDofMap(false)->getPOrder(elem, 0));

      const std::vector<Real>& JxW = fe->get_JxW();
      const std::vector<VectorFunction >& edge_phi = fe->getFunctions();
      const std::vector<Point>& xyz = fe->get_xyz();

      for (unsigned int qp=0; qp<qrule->n_points(); qp++) {
        const Point& point = xyz[qp];
        double Xden = 0.0;//std::abs(point(0)) < 2 ? 1 : 0;
        _exciton_sim->get_solution(elem, densityId, Xden, point);

        std::vector<double> Esolution;
        //std::map<ID, std::vector<double> >& solutions
        get_solution(elem, Efield + 0, Esolution, point);

	double epsilon;

	get_solution(elem, Epsilon, epsilon, point);
        Point E(Esolution[0], Esolution[1], Esolution[2]);

        E2 += JxW[qp] * (epsilon * E * E);

        if (inWell) {
          Fexc2[subdomain] += JxW[qp] * Xden;
          FexcE[subdomain] += JxW[qp] * std::sqrt(Xden) * E;
        }
      }
    }
  }

  delete fe;

  std::map<ID, double> Vi2;
  double sum_Vi = 0;

  //std::cout << "WC=" << Wc0 << " Wexc0=" << Wexc0 << " Wlt=" << Wlt << "\n";
  for (std::map<ID, double>::iterator it = Fexc2.begin(); it != Fexc2.end(); it++) {
    ID id = it->first;
    double a_b_GaN = 3e-9;
    double a_b_GaN_2D = 3e-9 / 2;

    double q = Wlt * a_b_GaN * a_b_GaN * a_b_GaN / (lwell * a_b_GaN_2D * a_b_GaN_2D);

    Vi2[id] = (Fexc2[id] == 0) ? 0.0 : (q * std::abs(Wc0) * ((FexcE[id] * FexcE[id]) / Fexc2[id] / E2));
    sum_Vi += Vi2[id];
  }

  Complex two(2, 0);


  WPolaritonLow = (Wc0 + Wexc0) / two - std::sqrt((Wc0 - Wexc0) / two * (Wc0 - Wexc0) / two + sum_Vi);
  std::cout << "Wpolariton low: " << WPolaritonLow << "\n";

  for (std::map<ID, double>::iterator it = Fexc2.begin(); it != Fexc2.end(); it++) {
    ID id = it->first;

    hopfieldCoeeficients[id] = Vi2[id] / (sum_Vi + std::abs((WPolaritonLow - Wexc0) * (WPolaritonLow - Wexc0)));
    //std::cout << "Hop " << id << " " << Vi2[id] << " " << sum_Vi << " " << std::abs((WPolaritonLow - Wexc0) * (WPolaritonLow - Wexc0)) << "\n";
    //std::cout << "Hop result: " << id << " " << hopfieldCoeeficients[id] << "\n";
  }
}
