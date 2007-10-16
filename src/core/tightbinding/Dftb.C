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
