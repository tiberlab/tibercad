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

Dftb::Dftb(void){
  _dftb_options.iPeriodic = 0;
  _dftb_options.nAtom = 0;
  _dftb_options.nType = 0;
};



Dftb::~Dftb(void){};


void 
Dftb::do_init(void){

  std::cerr << "Dftb Simulation Inizialization" << std::endl;


  // Getting reference to atomistic structure for calculation
  get_atomistic_structure();
  std::cerr << "Caught atomistic structure " << _atomistic_structure->get_name() << std::endl;

  // Building and searching SK files names
  build_names();

  build_structure_options();

  std::cerr << "Done " << std::endl;

}


void Dftb::do_solve(void){};

void Dftb::parse_options(void){};





void Dftb::build_names(void){

  std::string sk_name = "";
  std::ifstream file;
  std::vector<std::string> atom_types;
  atom_types = _atomistic_structure->get_atom_types ();



  // SK FILES NAMES

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
  _dftb_options.skNames[n_files * DFTBP_LC - 1] = '\0';


  //for (int i = 0; i < n_files * DFTBP_LC; i++) {std::cout << "Char " << i << " is " << _dftb_options.skNames[i] << std::endl;}
  //std::cout << "DFTB_OPTIONS.SKNAMES is " << _dftb_options.skNames << std::endl;




  // SPECIES NAMES
  _dftb_options.speciesNames = (char *) malloc( _atomistic_structure->N_types * DFTBP_MC * sizeof(char));

  counter = 0;
  for (int i = 0; i < _atomistic_structure->N_types; i++){

    for (int str_i = 0; str_i <  DFTBP_MC - 1; str_i++){
      if (str_i <  atom_types[i].size() ) _dftb_options.speciesNames[str_i + counter * DFTBP_MC] = atom_types[i][str_i];
      else _dftb_options.speciesNames[str_i +  counter * DFTBP_MC] = DFTBP_PADCHAR;
    }
    counter++;
  }

  _dftb_options.speciesNames[ _atomistic_structure->N_types * DFTBP_MC - 1] = '\0';


  //for (int i = 0; i <  _atomistic_structure->N_types * DFTBP_MC; i++) {std::cout << "Char " << i << " is " << _dftb_options.speciesNames[i] << std::endl;}


};



void Dftb::build_structure_options(){


  std::cout << "build_structure_options begin" << std::endl;

  _dftb_options.nAtom = _atomistic_structure->N_atoms;
  _dftb_options.nType = _atomistic_structure->N_types;

  _dftb_options.coords = new double[_dftb_options.nAtom * 3];

  std::vector<Atom> basis;
  basis = _atomistic_structure->get_structure_atoms();

  //! Setting coordinates in DFTB format
  for (int i = 0; i < _dftb_options.nAtom; i++){

    _dftb_options.coords[ (i*3) ] = basis[i].position(1);
    _dftb_options.coords[ (i*3) + 1 ] = basis[i].position(2);
    _dftb_options.coords[ (i*3) + 2 ] = basis[i].position(3);

  }

  //! Setting species in DFTB format
  _dftb_options.species = new int [ _dftb_options.nAtom ];

  for (int i = 0; i < _dftb_options.nAtom; i++){

    _dftb_options.species[i] = _atomistic_structure->get_type_index(basis[i].specie);

  }


 for (int i = 0; i < _dftb_options.nAtom * 3; i++){
   std::cout << "COORDS["<<i<<"] is " <<  _dftb_options.coords[ i ] << std::endl;
  }
 for (int i = 0; i < _dftb_options.nAtom; i++){
   std::cout << "SPECIE["<<i<<"] is " <<  _dftb_options.species[ i ] << std::endl;
  }


 _dftb_options.iPeriodic = _atomistic_structure->is_periodical;


};
