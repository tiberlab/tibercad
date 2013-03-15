// $Id: MaxwellAbsorption.C 2063 2010-09-03 13:11:49Z maufder $

#include "MaxwellAbsorption.h"
#include "SimulationEnvironment.h"
#include "Material.h"
#include "TiberMath.h"
#include "Database.h"
#include <assert.h>

#include "equation_systems.h"
#include "limits.h"
#include "Utils.h"
#include "point.h"
#include "math.h"
#include "Constants.h"
#include "SimulationInterface.h"
#include "Database.h"

#include "TiberModule.h"


MaxwellAbsorption::MaxwellAbsorption(const ModelOptions& options)
 : SimulationInterface(options)
{
//  has_solution_vector(false);
}


void MaxwellAbsorption::do_init() {
  std::string m = get_options().get_option("maxwell_simulation", "maxwell_boundary");
  maxwell = SimulationInterface::find_simulation(m);

  lambdaStart = get_options().get_option("lambdaStart", 400);
  lambdaEnd = get_options().get_option("lambdaEnd", 700);
  lambdaStep = get_options().get_option("lambdaStep", 10);

  plotAll = get_options().get_option("plotAll", true);
}

//=======================================================================================================//
void MaxwellAbsorption::do_solve() {
  std::string intencityFile = get_option("input_intencity",  Database::get_search_path() + "/optics/sun.dat");

  Database db2;
  db2.set_data_file(intencityFile);
  db2.get("lambda", input_lambda);
  db2.get("intensity", input_intencity);
  intensity_scaling = db2.get("scaling", 1.0);
  db2.close();

  int currentIndex = 0;

  const MeshBase& mesh = get_mesh();
  const MeshBase::const_element_iterator end_el =
                                  mesh.active_local_elements_end();

  ID id_efield = maxwell->get_solution_id("Efield");
  ID id_epsilon = maxwell->get_solution_id("Epsilon");
  ID id_epsilon_imag = maxwell->get_solution_id("Epsilon_imag");

  std::string outputPrefix = maxwell->get_name();

  for (double lambda = lambdaStart; lambda < lambdaEnd; lambda += lambdaStep) {
    while (currentIndex < input_lambda.size() && input_lambda[currentIndex] < lambda) {
      currentIndex++;
    }

    double currentIntensity = 1.0;
    if (lambda <= input_lambda[0]) {
      currentIntensity = input_intencity[0];
    } else if (lambda >= input_lambda[input_lambda.size() - 1]) {
      currentIntensity = input_intencity[input_lambda.size() - 1];
    } else {
      // linear interpolation

      currentIntensity = input_intencity[currentIndex-1] + (input_intencity[currentIndex] - input_intencity[currentIndex-1]) *
            (lambda - input_lambda[currentIndex - 1]) / (input_lambda[currentIndex] - input_lambda[currentIndex - 1]);

    }

    currentIntensity *= intensity_scaling;

    double W = 2 * M_PI / (lambda * 1e-9) * Constants::c;
    double W_in_eV = Constants::hbar * W / Constants::e;

    char buffer [50];
    sprintf(buffer, "%s_%.1lf", outputPrefix.c_str(), lambda);

    maxwell->get_options().set_option("W", W_in_eV);
    maxwell->get_options().set_option("output_prefix", std::string(buffer));
    maxwell->reinit();

    maxwell->solve();

    if (plotAll) {
      maxwell->plot();
    }

    double intencity = currentIntensity * lambdaStep;

    MeshBase::const_element_iterator el =
                                    mesh.active_local_elements_begin();

    // loop over all active elements
    for ( ; el != end_el ; ++el)
    {
      const Elem* elem = *el;

      bool firstTime = !solution_photon.count(elem->centroid());
      std::vector<Point> points;

      for (int i = 0; i < elem->n_nodes(); i++) {
        if (firstTime) {
          solution_energy[elem->centroid()].push_back(0);
          solution_photon[elem->centroid()].push_back(0);
        }
        points.push_back(elem->point(i));
      }

      std::map<ID, std::vector<double>> maxwell_solutions;
      // add sol vars to maxwell_solutions
      maxwell_solutions[0].resize(0);
      maxwell_solutions[1].resize(0);
      maxwell_solutions[2].resize(0);
      maxwell_solutions[3].resize(0);
      maxwell_solutions[4].resize(0);
      maxwell_solutions[5].resize(0);
      maxwell_solutions[6].resize(0);
      maxwell_solutions[7].resize(0);

      maxwell->get_solution(elem, maxwell_solutions, points, false);

      for (int i = 0; i < elem->n_nodes(); i++) {
        // We need absorbed count of photons per m^3:     I * w/c * Imag(epsilon) / hw
        // It should be multiplied by input intencity ofc.

        std::complex<double> n(maxwell_solutions[id_epsilon][i], maxwell_solutions[id_epsilon_imag][i]);
        n = std::sqrt(n);

        double eSquared = maxwell_solutions[id_efield][3*i] * maxwell_solutions[id_efield][3*i] +
            maxwell_solutions[id_efield][3*i+1] * maxwell_solutions[id_efield][3*i+1] +
            maxwell_solutions[id_efield][3*i+2] * maxwell_solutions[id_efield][3*i+2];

        double realESquared = eSquared * intencity / (Constants::epsilon * Constants::c / 2);

        //double absrobed_I = 2 * realESquared * 2 * n.imag() * W /Constants::c / 100; //
        double absrobed_I = n.imag() * Constants::epsilon * realESquared * n.real() * W / 1e6; // 1e6 m-3 -> cm-3
        double absorbed_ph = absrobed_I / Constants::hbar / W;

        solution_energy[elem->centroid()][i] += absrobed_I;
        solution_photon[elem->centroid()][i] += absorbed_ph;
      }
    }
  }

}


void
MaxwellAbsorption::do_setup_solution_variables(void)
{
  declare_solution_ext("Absorbed energy", Absorption_Energy, SolutionDescriptor::REAL, SolutionDescriptor::NODES, "J/cm^3/c");
  declare_solution_ext("Absorbed photons", Absorption_Photon, SolutionDescriptor::REAL, SolutionDescriptor::NODES, "1/cm^3/c");
}

void MaxwellAbsorption::get_solution_secure(const Elem* elem,
    std::map<ID, std::vector<double> >& solutions,
    const std::vector<Point>& points)
{
  if (solutions.count(Absorption_Energy) || solutions.count(Absorption_Photon)) {

    std::vector<double>& sol_e = solutions[Absorption_Energy];
    sol_e.resize(points.size());

    std::vector<double>& sol_ph = solutions[Absorption_Photon];
    sol_ph.resize(points.size());

    for (int i = 0; i < points.size(); i++) {
      sol_e[i] = 0;
      sol_ph[i] = 0;
    }

    const unsigned int dim = get_mesh().mesh_dimension();

    AutoPtr<FEBase> fe(build_finite_element(dim, FEType()));

    // element shape functions
    const std::vector<std::vector<double> >& phi = fe->get_phi();

    fe->reinit(elem, &points);

    for (int i = 0; i < elem->n_nodes(); i++) {
      for (int qp = 0; qp < points.size(); qp++) {
        sol_e[qp] += phi[i][qp] * solution_energy[elem->centroid()][i];
        sol_ph[qp] += phi[i][qp] * solution_photon[elem->centroid()][i];
      }
    }
  }
}
