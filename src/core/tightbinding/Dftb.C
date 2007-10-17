//modules includes
#include "Dftb.h"

#include "BoundaryProperties.h"
#include "PhysicalModel.h"
#include "SimulationOptions.h"
#include "TightBindingModel.h"
#include "ElectricalContact.h"
#include "SimulationOptions.h"

using namespace std;

//--------------------------------------------------------------

Dftb::Dftb(void){};

Dftb::~Dftb(void){};

void Dftb::do_init(void){

  cerr << "Dftb Simulation Inizialization" << endl;

  _atomistic_structure = NULL;

  if (get_options().find_option("atomistic_structure") )
    {
      string name;
      name = get_options().get_option("atomistic_structure", "none");
      if (name.compare("none") != 0){
      _atomistic_structure = get_environment().get_device().get_atomistic_structure(name);
      }
    }
  else 
    {
     cerr << "ERROR in Tight Binding Simulation: an atomistic structure  must be specified " 
      << get_name() << endl;
     exit(0);
    }


}


void Dftb::do_solve(void){};

void Dftb::parse_options(void){};
