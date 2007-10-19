//modules includes
#include "Dftb.h"
#include "BoundaryProperties.h"
#include "PhysicalModel.h"
#include "SimulationOptions.h"
#include "TightBindingModel.h"
#include "ElectricalContact.h"
#include "SimulationOptions.h"
#include "DftbpWrapper.h"
#include "dftbp.h"

#include <fstream>


//--------------------------------------------------------------

Dftb::Dftb(void){};



Dftb::~Dftb(void){};


void 
Dftb::do_init(void){

  std::cerr << "Dftb Simulation Inizialization" << std::endl;


  // Getting reference to atomistic structure for calculation
  get_atomistic_structure();
  std::cerr << "Caught atomistic structure " << _atomistic_structure->get_name() << std::endl;

  // Building and searching SK files names
  build_sk_names();

  std::cerr << "Done " << std::endl;

}


void Dftb::do_solve(void){};

void Dftb::parse_options(void){};



void Dftb::build_sk_names(void){

  std::string sk_name = "";
  std::ifstream file;
  std::vector<std::string> atom_types;
  atom_types = _atomistic_structure->get_atom_types ();

  //! This is the SK files path
  const std::string prefix = "";

  // Iterative factorial is enough, number of species is a small value
  int n_files = 0;
  n_files = _atomistic_structure->N_types * _atomistic_structure->N_types;
 
  std::cout << "NUMBER OF FILES IS " << n_files << std::endl;

  //Static allocation. Size must be decided according to dftbp.h parameters
  _dftb_options.skNames = (char *) malloc(n_files * DFTBP_LC * sizeof(char));

  // Putting NULL character at the end!!!!!!!!!!!!!! Don't know if really needed by DFTB interface!!!!!!!!!!!!!!!!!!
  //memset( _dftb_options.skNames, '\0', n_files * DFTBP_LC );


  // Cycle upon species and build names
  int counter = 0;
  for (int i = 0; i < ( _atomistic_structure->N_types); i++){
    for (int j = 0; j < (_atomistic_structure->N_types); j++){
    
      sk_name.clear();
      sk_name.append(prefix);
      sk_name.append(atom_types[i]); sk_name.append("-"); 
      sk_name.append(atom_types[j]); sk_name.append(".skf");
      std::cout << "Checking sk_name " << sk_name << std::endl;

      if (i == j){
	file.open(sk_name.c_str());
	if ( !(file.is_open()) ) {std::cerr << "ERROR IN DFTB: COULD NOT FIND SK FILE " 
					    << sk_name << std::endl;}
	file.close();
      }

      if (i != j){

	file.open(sk_name.c_str());

	if ( !(file.is_open()) ) {
	  file.close();
	  sk_name.clear();
	  sk_name.append(prefix);
	  sk_name.append(atom_types[j]); sk_name.append("-"); 
	  sk_name.append(atom_types[i]); sk_name.append(".skf");
	}
	file.close();

	file.open(sk_name.c_str());
	if ( !(file.is_open()) ) {std::cerr << "ERROR IN DFTB: COULD NOT FIND SK FILE " 
					    << sk_name << std::endl;}

	file.close();

      }

      if ( sk_name.size() > DFTBP_LC ) {std::cerr << "ERROR IN DFTB: SK FILENAME " << sk_name 
						  <<" IS TOO LONG " << std::endl;}

      std::cout << "SKNAME IS " << sk_name << std::endl;
      std::cout << "Dimension is " << sk_name.size() << std::endl;
      std::cout << "Counter is " << counter << std::endl;

      for (int str_i = 0; str_i <  DFTBP_LC - 1; str_i++){
	if (str_i < sk_name.size()) _dftb_options.skNames[str_i + counter * DFTBP_LC] = sk_name[str_i];
	else _dftb_options.skNames[str_i +  counter * DFTBP_LC] = DFTBP_PADCHAR;
      }

  counter++;
      
    }
  }
  _dftb_options.skNames[n_files * DFTBP_LC - 1] = NULL;

  //for (int i = 0; i < n_files * DFTBP_LC; i++) {std::cout << "Char " << i << " is " << _dftb_options.skNames[i] << std::endl;}
  std::cout << "DFTB_OPTIONS.SKNAMES is " << _dftb_options.skNames << std::endl;

};
