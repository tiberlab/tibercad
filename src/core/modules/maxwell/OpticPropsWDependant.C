
#include "OpticPropsWDependant.h"

#include "SimulationInterface.h"
#include "Material.h"


TIBER_MODULE(OpticPropsWDependant, opticmodel, w_dependant)

void OpticPropsWDependant::read_database() {
  const Database& db = get_database();
  db.set_section("permittivity");

  dataFile = db.get("dataFile", "");
}

void OpticPropsWDependant::do_init(void) {
  get_parameter("dataFile", dataFile);

  if (dataFile == "") {
    throw InitFailedException("Empty file name");
  } else {
    const std::string path = Database::get_search_path() + "/optics/" + get_database().get_material() + ".dat";
    Database db2(get_database());
    db2.set_data_file(path);

    db2.get("W", params);
    db2.get("eps", epsilon);
    db2.get("eps_imag", epsilon_imag);

  }
}

void OpticPropsWDependant::set_param(double param) {
  int newIndex = currentIndex;

  while (params[newIndex] < param && newIndex < params.size()) {
    newIndex++;
  }

  currentIndex = (newIndex < params.size()) ? newIndex : params.size() - 1;
}
