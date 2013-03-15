
#include "OpticPropsWDependant.h"

#include "SimulationInterface.h"
#include "Material.h"
#include "Messages.h"
#include "Constants.h"
#include "MaxwellBoundaryEquations.h"

#include "TiberModule.h"



void OpticPropsWDependant::read_database() {

  const Database& db = get_database();
  db.set_section("permittivity");

  dataFile = db.get("optic_datafile", Database::get_search_path() + "/optics/" + db.get_material() + ".dat");

}

void OpticPropsWDependant::do_init(void) {

  currentIndex = 0;

  get_parameter("optic_datafile", dataFile);

  if (dataFile == "") {
    throw InitFailedException("Empty file name");
  } else {

    Database db2;
    db2.set_data_file(dataFile);


    db2.get("lambda", params);

    std::vector<double> epsilon_real;
    std::vector<double> epsilon_imag;

    db2.get("eps", epsilon_real);
    db2.get("eps_imag", epsilon_imag);

    bool useEpsilon = epsilon_real.size() > 0;

    if (!useEpsilon) {
    	db2.get("n", epsilon_real);
	db2.get("k", epsilon_imag);
    }

    epsilon.resize(epsilon_real.size());

    for (int i = 0; i < epsilon.size(); i++) {
      if (useEpsilon) {
      	epsilon[i] = Complex(epsilon_real[i], epsilon_imag[i]);
      } else {
        epsilon[i] = Complex(epsilon_real[i], epsilon_imag[i]) * Complex(epsilon_real[i], epsilon_imag[i]);
      }
    }

    db2.close();

  }
}

void OpticPropsWDependant::do_reinit() {

	//std::cout << "reinit\n";
  MaxwellBoundaryEquations* maxwell = dynamic_cast<MaxwellBoundaryEquations*>(SimulationInterface::get_simulation(get_simulator_id()));

  double W = maxwell->getW();

  double lambda =  2 * M_PI * Constants::c / W / 1.0e-9; // in nm

  int newIndex = currentIndex;

  while (newIndex < params.size() && params[newIndex] < lambda) {
    newIndex++;
  }

  currentIndex = (newIndex < params.size()) ? newIndex : (params.size() - 1);
  
  if (lambda <= params[0]) {
    currentEpsilon = epsilon[0];
  } else if (lambda >= params[params.size() - 1]) {
    currentEpsilon = epsilon[epsilon.size() - 1];
  } else {
    // linear interpolation
    currentEpsilon = epsilon[currentIndex-1] + (epsilon[currentIndex] - epsilon[currentIndex-1]) *
          (lambda - params[currentIndex - 1]) / (params[currentIndex] - params[currentIndex - 1]);
  }
}

Complex OpticPropsWDependant::get_dielectric_constant() const {
  return currentEpsilon;
}
