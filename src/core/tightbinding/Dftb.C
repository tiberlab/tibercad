// $Id$

#include "tiber_config.h"

#ifdef ENABLE_DFTB

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

Dftb::Dftb(void)
  : _dftb_options()
{
  inst = new DftbpWrapper;
};


Dftb::~Dftb(void){};


Dftb::DftbOptions::DftbOptions(void){
  iPeriodic = 0;
  nAtom = 0;
  nType = 0;
  latVecs = new double[9];
  for (int i = 0; i < 9; i++) {latVecs[i] = 0.0;}
  samplingcoeffs = new double[9];
  samplingshift = new double[3];
};


Dftb::DftbOptions::~DftbOptions(void){
  delete samplingcoeffs;
  delete samplingshift;
};


void Dftb::parse_options(void){

  //Set up file names for Slater Koster files
  build_names();

  //Set up structure options from atomistic structure
  build_structure_options();

  //Set up input options for building up simulation
  build_input_options();

  //Set up input options for the solver (up to now only internal solver is supported)
  _dftb_solver_options.solver = "external";


#ifdef DEBUG
  print_dftb_options();
#endif

};


void 
Dftb::do_init(void){

  std::cerr << "Dftb Simulation Inizialization" << std::endl;


  // Getting reference to atomistic structure for calculation
  get_atomistic_structure();
  std::cerr << "Caught atomistic structure " << _atomistic_structure->get_name() << std::endl;


  // Setting options for DFTB+ calls
  parse_options();


  //  Initialize Dftb instance it and sets parameters
  std::cout << "Am I periodic?? " << _dftb_options.iPeriodic << std::endl;
  inst->fill_param(_dftb_options.nAtom, _dftb_options.nType, 
		   _dftb_options.eTemp, _dftb_options.iPeriodic, _dftb_options.speciesNames, 
		   _dftb_options.species);

  std::cout << "fill parameter done" << std::endl;

  inst->addskdata(_dftb_options.skNames, _dftb_options.mAngs, 
		  _dftb_options.orbResolved, _dftb_options.skInterp, _dftb_options.nType);

  std::cout << "addskdata done" << std::endl;

  if (_dftb_options.iPeriodic == 1) 
    {inst->addlattice(_dftb_options.latVecs);
      std::cout << "addlattice done" << std::endl;}

  if (_dftb_options.iPeriodic == 1) 
    {if (_dftb_options.supersampling = true) 
	{
	  //Parameter noinv is set to 1 (use inversion simmetry). 
	  //To 0 will generate the complete set (both positive and negative k points)
	  int noinv = 1;
	  inst->addsupersampling(_dftb_options.samplingcoeffs, _dftb_options.samplingshift, noinv);
	}
else
      inst->addkpoints(_dftb_options.nkPoints, _dftb_options.kPoints, _dftb_options.kWeights);
      std::cout << "addkPoints done" << std::endl;}

  std::cout << "initdftb begins" << std::endl;
 
  inst->initdftb();

  std::cout << "initdftb done" << std::endl; 

  inst->up_coords(_dftb_options.nAtom, _dftb_options.coords);

  std::cout << "up_coords" << std::endl; 
 
 
}



void Dftb::do_solve(void){

#ifdef DEBUG
  std::cout << "Calling Dftb->do_solve() " << std::endl;
#endif

  //get_energy is anyway needed for the correct assembling of Hamiltonian
   double energy = 0;
   inst->get_energy(energy);

  if (_dftb_solver_options.solver.compare("internal") == 0) {

   double* charges;
  charges = new double[_dftb_options.nAtom];
  inst->getchargesperatom(_dftb_options.nAtom, charges);
 
  _atomistic_structure->print_structure("TB_out.xyz",charges);

  }

  int nrow,ncol,nzval,isreal;
  //inst->getrealhamiltonian(nrow,ncol,nzval,isreal);

  int *colind, *rowpnt;
  double* val;
  std::string matrix;
  matrix = "H";
  double* kPoint;
  kPoint = new double[3];
  kPoint[0] = 0.0; kPoint[1] = 0.0; kPoint[2] = 0.0;
  inst->getmatrix(nrow, ncol, nzval, isreal, colind, rowpnt, val, matrix);

#ifdef DEBUG
  std::cout << "Dfbt->solve() done" << std::endl;
#endif

};







void Dftb::build_names(void){

  std::string sk_name = "";
  std::ifstream file;
  std::vector<std::string> atom_types(_atomistic_structure->get_atom_types());



  // SK FILES NAMES

  //! This is the SK files path
  const std::string prefix = "";

  int n_files = 0;
  n_files = _atomistic_structure->get_N_types() * _atomistic_structure->get_N_types();


  //Static allocation. Size must be decided according to dftbp.h parameters
  _dftb_options.skNames = (char *) malloc(n_files * DFTBP_LC * sizeof(char));

  // Putting NULL character at the end!!!!!!!!!!!!!! Don't know if really needed by DFTB interface!!!!!!!!!!!!!!!!!!
  //memset( _dftb_options.skNames, '\0', n_files * DFTBP_LC );


  // Cycle upon species and build names
  int counter = 0;
  for (int i = 0; i < ( _atomistic_structure->get_N_types() ); i++){
    for (int j = 0; j < (_atomistic_structure->get_N_types() ); j++){
    
      sk_name.clear();
      sk_name.append(prefix);
      sk_name.append(atom_types[i]); sk_name.append("-"); 
      sk_name.append(atom_types[j]); sk_name.append(".skf");

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

      for (int str_i = 0; str_i <  DFTBP_LC - 1; str_i++){
	if (str_i < sk_name.size()) _dftb_options.skNames[str_i + counter * DFTBP_LC] = sk_name[str_i];
	else _dftb_options.skNames[str_i +  counter * DFTBP_LC] = DFTBP_PADCHAR;
      }
      counter++;
    }
  }
  _dftb_options.skNames[n_files * DFTBP_LC - 1] = '\0';


  // SPECIES NAMES
  _dftb_options.speciesNames = (char *) malloc( _atomistic_structure->get_N_types() * DFTBP_MC * sizeof(char));

  counter = 0;
  for (int i = 0; i < _atomistic_structure->get_N_types(); i++){

    for (int str_i = 0; str_i <  DFTBP_MC - 1; str_i++){
      if (str_i <  atom_types[i].size() ) _dftb_options.speciesNames[str_i + counter * DFTBP_MC] = atom_types[i][str_i];
      else _dftb_options.speciesNames[str_i +  counter * DFTBP_MC] = DFTBP_PADCHAR;
    }
    counter++;
  }

  _dftb_options.speciesNames[ _atomistic_structure->get_N_types() * DFTBP_MC - 1] = '\0';


};



void Dftb::build_structure_options(){


  std::cout << "build_structure_options begin" << std::endl;

  _dftb_options.nAtom = _atomistic_structure->get_N_atoms();
  _dftb_options.nType = _atomistic_structure->get_N_types();

  _dftb_options.coords = new double[_dftb_options.nAtom * 3];

  std::vector<Atom> basis;
  basis = _atomistic_structure->get_structure_atoms();

  //! Setting coordinates in DFTB format
  for (int i = 0; i < _dftb_options.nAtom; i++){

    _dftb_options.coords[ (i*3) ] = basis[i].get_position(1);
    _dftb_options.coords[ (i*3) + 1 ] = basis[i].get_position(2);
    _dftb_options.coords[ (i*3) + 2 ] = basis[i].get_position(3);

  }

  //! Setting species in DFTB format
  _dftb_options.species = new int [ _dftb_options.nAtom ];

  for (int i = 0; i < _dftb_options.nAtom; i++){
    _dftb_options.species[i] = _atomistic_structure->get_type_index(basis[i].get_specie());

  }

  _dftb_options.iPeriodic = _atomistic_structure->is_periodic();

  //If system is not periodical lattice vectors are set to 0 by default in AtomisticStructure.
  //Check if huge value is needed instead of zero value
  double* tmp_periodicity_vectors;
  _dftb_options.latVecs = _atomistic_structure->get_periodicity_vectors();

  //If lattice vectors are found in input file (in order x1, y1, z1, x2...), overwrite them
  if ( get_options().find_option("latvecs") ) {

    //If lattice vectors are specified, structure is considered periodical
    _dftb_options.iPeriodic = 1;
    std::vector<double> tmp_latvecs;
    get_options().get_option("latvecs", tmp_latvecs);
    if (tmp_latvecs.size() != 9) std::cerr << "ERROR: latvecs must be a 9 doubles array" << std::endl;
    for (int i = 0; i < 9; i++) {_dftb_options.latVecs[i] = tmp_latvecs[i];}

  }


};


void Dftb::build_input_options(){

  _dftb_options.eTemp = get_options().get_option("eTemp", 1e-8);

  _dftb_options.orbResolved = get_options().get_option("orbResolved", 0);

  _dftb_options.skInterp = get_options().get_option("skInterp", 2);

  if (get_options().find_option("supercellfolding")) {
    _dftb_options.iPeriodic = 1;
    _dftb_options.supersampling = true;
    std::vector<double> supersamplingdata;

    get_options().get_option("supercellfolding", supersamplingdata);

    if (supersamplingdata.size() != 12) std::cerr << "Warning: supercellfolding options must be 12 numbers (9 for supercell vecotrs coefficients and 3 for shift)" << std::endl;

    for (int i = 0; i < 9 ; i++) {
      _dftb_options.samplingcoeffs[i] = supersamplingdata[i];
    }

 for (int i = 9; i < 12 ; i++) {
      _dftb_options.samplingshift[i - 9] = supersamplingdata[i];
    }
    
  }


  if (get_options().find_option("iPeriodic")){
    _dftb_options.iPeriodic = get_options().get_option("iPeriodic", 0);
  }



  //CHANGE THE METHOD TO ASSIGN MAXIMUM ANGOLAR MOMENTUM!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  _dftb_options.mAngs = new int[_dftb_options.nType];
  for (int i = 0; i < _dftb_options.nType; i++)   _dftb_options.mAngs[i] = 1;


  if (_dftb_options.iPeriodic == 1){
    //CAREFULL!!! TEMPORARY DEFAULTS! RIGHT ONES MUST BE DISCUSSED
    _dftb_options.nkPoints = get_options().get_option("nkPoints", 4);


  


    //A METHOD FOR INSERTING KPOINTS MUST STILL BE DECIDED. ONLY DEFAULT AVAILABLE
    _dftb_options.kPoints = new double [   _dftb_options.nkPoints * 3 ];
    if ( _dftb_options.nkPoints == 1) {
      _dftb_options.kPoints[0] = 0.25;
      _dftb_options.kPoints[1] = 0.25;
      _dftb_options.kPoints[2] = 0.25;
    }
    else
      if ( _dftb_options.nkPoints == 4) {
	_dftb_options.kPoints[0] = 0.25;
	_dftb_options.kPoints[1] = 0.25;
	_dftb_options.kPoints[2] = 0.25;
	_dftb_options.kPoints[3] = -0.25;
	_dftb_options.kPoints[4] = 0.25;
	_dftb_options.kPoints[5] = 0.25;
	_dftb_options.kPoints[6] = 0.25;
	_dftb_options.kPoints[7] = -0.25;
	_dftb_options.kPoints[8] = 0.25;
	_dftb_options.kPoints[9] = -0.25;
	_dftb_options.kPoints[10] = -0.25;
	_dftb_options.kPoints[11] = 0.25;
 
      }
      else 
	if ( _dftb_options.nkPoints == 8) {
	  _dftb_options.kPoints[0] = 0.25;
	  _dftb_options.kPoints[1] = 0.25;
	  _dftb_options.kPoints[2] = 0.25;
	  _dftb_options.kPoints[3] = -0.25;
	  _dftb_options.kPoints[4] = 0.25;
	  _dftb_options.kPoints[5] = 0.25;
	  _dftb_options.kPoints[6] = 0.25;
	  _dftb_options.kPoints[7] = -0.25;
	  _dftb_options.kPoints[8] = 0.25;
	  _dftb_options.kPoints[9] = -0.25;
	  _dftb_options.kPoints[10] = -0.25;
	  _dftb_options.kPoints[11] = 0.25;
	  _dftb_options.kPoints[12] = 0.5;
	  _dftb_options.kPoints[13] = 0.5;
	  _dftb_options.kPoints[14] = 0.5;
	  _dftb_options.kPoints[15] = -0.5;
	  _dftb_options.kPoints[16] = 0.5;
	  _dftb_options.kPoints[17] = 0.5;
	  _dftb_options.kPoints[18] = 0.5;
	  _dftb_options.kPoints[19] = -0.5;
	  _dftb_options.kPoints[20] = 0.5;
	  _dftb_options.kPoints[21] = -0.5;
	  _dftb_options.kPoints[22] = -0.5;
	  _dftb_options.kPoints[23] = 0.5;
	}
	else {std::cout << "Error, up to now only default 4 k points available" << std::endl;exit(0);}
  
    _dftb_options.kWeights = new double[  _dftb_options.nkPoints];
    for (int i = 0; i < _dftb_options.nkPoints; i++){ 
      _dftb_options.kWeights[i] = 1.0;
    }

  }

};


void Dftb::print_dftb_options(void){

  std::cout << "DFTB_OPTIONS: " << std::endl;

  // std::cout << "coords are " << std::endl;
//   for (int i = 0; i < _dftb_options.nAtom * 3; i++){
//     std::cout << "coords["<<i<<"] is " <<  _dftb_options.coords[ i ] << std::endl;
//   }

//   std::cout << "species are " << std::endl;
//   for (int i = 0; i < _dftb_options.nAtom; i++){
//     std::cout << "species["<<i<<"] is " <<  _dftb_options.species[ i ] << std::endl;
//   }

  int n_files = 0;
  // n_files = _atomistic_structure->N_types * _atomistic_structure->N_types;
 //  std::cout << "skNames are " << std::endl;
//   for (int i = 0; i < n_files * DFTBP_LC; i++) {std::cout << "Char " << i << " is " << _dftb_options.skNames[i] << std::endl;}
//   std::cout << "skNames string lenght is " << strlen(_dftb_options.skNames);
 
//   std::cout << "speciesNames are " << std::endl;
//   for (int i = 0; i <  _atomistic_structure->N_types * DFTBP_MC; i++) {std::cout << "Char " << i << " is " << _dftb_options.speciesNames[i] << std::endl;}
//   std::cout << "speciesNames string lenght is " << strlen(_dftb_options.speciesNames);

  std::cout << "latVecs are " << std::endl;
  for (int i = 0; i < 9; i++){
    std::cout << "latvecs[" <<i<<"]"<< _dftb_options.latVecs[i] << std::endl; 
  }

  std::cout << "iPeriodic is " << _dftb_options.iPeriodic << std::endl;
  std::cout << "eTemp is " << _dftb_options.eTemp << std::endl;
  std::cout << "nAtom is " << _dftb_options.nAtom << std::endl;
  std::cout << "nType is " << _dftb_options.nType << std::endl;
  std::cout << "skInterp is " << _dftb_options.skInterp << std::endl;

};


#endif // ENABLE_DFTB
