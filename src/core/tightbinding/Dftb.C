//modules includes
#include "Dftb.h"
#include "BoundaryProperties.h"
#include "PhysicalModel.h"
#include "SimulationOptions.h"
#include "TightBindingModel.h"
#include "ElectricalContact.h"
#include "SimulationOptions.h"

#include <cmath>

//--------------------------------------------------------------

Dftb::Dftb(void){
  _atomistic_structure = NULL;
};



Dftb::~Dftb(void){};


void 
Dftb::do_init(void){

  std::cerr << "Dftb Simulation Inizialization" << std::endl;


  // Getting reference to atomistic structure for calculation
  _atomistic_structure = get_atomistic_structure();

  std::cerr << "Caught atomistic structure " << _atomistic_structure->get_name() << std::endl;

  // Building and searching SK files names
  char* a = " ";
  a = build_sk_names();


}


void Dftb::do_solve(void){};

void Dftb::parse_options(void){};


char* Dftb::build_sk_names(void){

  std::string sk_name;

  //! This is the SK files path
  const std::string path = " ";

  // Iterative factorial is enough, number of species is a small value
  int n_files = 0;
  n_files = _atomistic_structure->N_types * _atomistic_structure->N_types;
 
  std::cout << "NUMBER OF FILES IS " << n_files << std::endl;

  //Static allocation. Size must be decided according to dftbp.h parameters



};
