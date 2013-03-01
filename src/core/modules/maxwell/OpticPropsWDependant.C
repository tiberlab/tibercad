
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

  dataFile = db.get("optic_datafile", "optics/" + db.get_material() + ".dat");

}

void OpticPropsWDependant::do_init(void) {

  currentIndex = 0;

  get_parameter("optic_datafile", dataFile);

  if (dataFile == "") {
    throw InitFailedException("Empty file name");
  } else {


    const std::string path = Database::get_search_path() + "/" + dataFile;
    Database db2(get_database());
    db2.set_data_file(path);

    db2.get("lambda", params);
    db2.get("eps", epsilon);
    db2.get("eps_imag", epsilon_imag);

    db2.close();

  }
}

void OpticPropsWDependant::do_reinit() {
  MaxwellBoundaryEquations* maxwell = dynamic_cast<MaxwellBoundaryEquations*>(SimulationInterface::get_simulation(get_simulator_id()));

  double W = maxwell->getW();

  double lambda =  2 * M_PI * Constants::c / W;

  int newIndex = currentIndex;

  while (newIndex < params.size() && params[newIndex] < lambda) {
    newIndex++;
  }

  currentIndex = (newIndex < params.size()) ? newIndex : params.size() - 1;
}
