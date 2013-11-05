#include "tiber_config.h"

#ifdef ENABLE_DFTB

//modules includes
#include "Dftb.h"
#include "BoundaryProperties.h"
#include "PhysicalModel.h"
#include "SimulationOptions.h"
#include "DftbModel.h"
#include "SimulationOptions.h"
#include "DftbpWrapper.h"
#include "dftbp.h"
#include "mesh.h"

#include <fstream>
#include <sstream>
#include <utility>
using namespace std;


//--------------------------------------------------------------

Dftb::Dftb(const ModelOptions& options)
: TightBinding(options),
  _dftb_options()
{
  inst = new DftbpWrapper;
  _max_shell = 0;
};


Dftb::~Dftb(void)
{
  delete inst;
  inst = NULL;
};


Dftb::DftbOptions::DftbOptions(void)
:iPeriodic(0),
nAtom(0),
nType(0),
supersampling(false)
{
  latVecs = new double[9];
  for (int i = 0; i < 9; i++) {latVecs[i] = 0.0;}
  samplingcoeffs = new double[9];
  samplingshift = new double[3];
};


Dftb::DftbOptions::~DftbOptions(void)
{
  //TODO: perche' crasha con questa roba nel codice scommentata???
  //delete[] samplingcoeffs;
  //delete[] samplingshift;
  //delete[] latVecs;
};


void Dftb::parse_options(void){

  std::ifstream file;
  std::string specie;
  int tmp_shell;
  std::string line, record;

  //====================================
  //Specifies where to search for Sk files
  //TODO: the right set should be provided as a physical model!
  //TODO: uncomment when a SK database exists in tibercad and add to string prefix
  //const std::string database_path = Database::get_default_search_path();
  //const std::string database_path = Database::get_default_search_path();

  const std::string sk_path = "./SK/infos.dat";

  //Informations about parametrization are expected to be inside a file called infos.dat
  //with the following style:
  // <specie>   <maximum angular momentum + 1>
  // # is considered to start a comment line


  file.open(sk_path.c_str(), std::ifstream::in);

  if (!file)
    {
      std::cerr << "Unable to open file " << sk_path << ". Cannot load Slater Koster parametrization infos in DFTB \n";
      exit(1);   // call system to stop
    }

  while (getline(file, line))
    {
      std::stringstream line_string(line);

      line_string >> record;

      if (record.substr(0,1).compare("#") == 0) continue;

      line_string >> tmp_shell;
      _shell[record] = tmp_shell;

      if (_max_shell < tmp_shell) _max_shell = tmp_shell;

    }

  //=====================================




  //Set up file names for Slater Koster files
  build_names();

  //Set up structure options from atomistic structure
  build_structure_options();

  //Set up input options for building up simulation
  build_input_options();

  //Set up input options for the solver (up to now only internal solver is supported)
  _dftb_solver_options.solver = "internal";


#ifdef DEBUG
  print_dftb_options();
#endif

};


void
Dftb::do_init(void){

  std::cerr << "Dftb Simulation Initialisation" << std::endl;

  TightBinding::do_init();

  // Setting options for DFTB+ calls
  parse_options();
  //  Initialize Dftb instance it and sets parameters
  inst->fill_param(_dftb_options.nAtom, _dftb_options.nType,
      _dftb_options.eTemp, _dftb_options.iPeriodic, _dftb_options.specieNames,
      _dftb_options.species);


  inst->addskdata(_dftb_options.skNames, _dftb_options.mAngs,
      _dftb_options.orbResolved, _dftb_options.skInterp, _dftb_options.nType);


  obtain_hubbard_parameters();


  if (_dftb_options.iPeriodic == 1)
    {
      inst->addlattice(_dftb_options.latVecs);

      if (_dftb_options.iPeriodic == 1)
        {
          if (_dftb_options.supersampling == true)
            {
              //Parameter noinv is set by default to 1 (use inversion simmetry).
              //To 0 will generate the complete set (both positive and negative k points)
              int noinv;
              noinv = get_options().get_option("noinv", 1);
              inst->addsupersampling(_dftb_options.samplingcoeffs, _dftb_options.samplingshift, noinv);
            }
          else
            {
              inst->addkpoints(_dftb_options.nkPoints, _dftb_options.kPoints, _dftb_options.kWeights);
            }
        }
    }

  std::cout << "initdftb begins" << std::endl;

  inst->initdftb();

  std::cout << "initdftb done" << std::endl;

  inst->up_coords(_dftb_options.nAtom, _dftb_options.coords);

  //Check for external potential

  if (get_options().find_option("potential_simulation"))
    {
      _dftbp_coupling_options.potential_sim = get_options().get_option("potential_simulation", "no_sim");
      add_shifts();
    }
}


void Dftb::do_solve(void){

#ifdef DEBUG
  std::cout << "Calling Dftb->do_solve() " << std::endl;
#endif

  print_dftb_options();

  // N.B. get_energy is anyway needed for the correct assembling of Hamiltonian,
  // but for the SCC module it has to be in do_solve
  double energy = 0;
  inst->get_energy(energy);

  std::cout << "DFTB library computed a total system energy of: " << energy << std::endl;

  //TODO: non so come funziona il settaggio delle grandezze da calcolare quando
  //queste siano richieste da altri e non siano prettamente dati di output

  _dftb_solver_options.solver="internal";
  if (_dftb_solver_options.solver.compare("internal") == 0) {

    double* charges;
    charges = new double[_dftb_options.nAtom];
    inst->getnetchargesperatom(_dftb_options.nAtom, charges);

    _mulliken_netcharges.resize(_dftb_options.nAtom);
    for (unsigned int i = 0; i <_dftb_options.nAtom; i++ ) _mulliken_netcharges[i] = charges[i];

    for (unsigned int i=0; i<_mulliken_netcharges.size(); i++) std::cout << _mulliken_netcharges[i];
    _atomistic_structure->print_structure("TB_out.xyz",charges);

  }

  //std::cout << "charge in middle point is " << build_rho(1.3, 1.3, 1.3) << std::endl;
  //std::cout << "charge in 0 is " << build_rho(0.0, 0.0, 0.0) << std::endl;
  //int nrow,ncol,nzval,isreal;
  //inst->getrealhamiltonian(nrow,ncol,nzval,isreal);

  //int *colind, *rowpnt;
  //double* val;
  //std::string matrix;
  //matrix = "H";
  //double* kPoint;
  //inst->getmatrix(nrow, ncol, nzval, isreal, colind, rowpnt, val, matrix);


#ifdef DEBUG
  std::cout << "Dfbt->solve() done" << std::endl;
#endif

};







void Dftb::build_names(void){

  std::string sk_name = "";
  std::ifstream file;

  _dftb_options.specieNameStrings.clear();
  _dftb_options.specieNameStrings.resize( _atomistic_structure->get_atom_types().size() );

  for (unsigned int i = 0; i < _atomistic_structure->get_atom_types().size(); i++)
    {
      _dftb_options.specieNameStrings[i] = (_atomistic_structure->get_atom_types()[i]);
      std::cout << "specie name is " << _dftb_options.specieNameStrings[i] << std::endl;
    }


  // SK FILES NAMES

  //! This is the SK files path

  //TODO: uncomment when a SK database exists in tibercad and add to string prefix
  //const std::string database_path = Database::get_default_search_path();
  const std::string prefix = "./SK/";


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
      sk_name.append( _dftb_options.specieNameStrings[i]); sk_name.append("-");
      sk_name.append( _dftb_options.specieNameStrings[j]); sk_name.append(".skf");

      //TODO:default local directory is used: to change
      //sk_name = "SK/" + sk_name;

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
          sk_name.append( _dftb_options.specieNameStrings[j]); sk_name.append("-");
          sk_name.append( _dftb_options.specieNameStrings[i]); sk_name.append(".skf");
        }
        file.close();

        file.open(sk_name.c_str());
        if ( !(file.is_open()) ) {std::cerr << "ERROR IN DFTB: COULD NOT FIND SK FILE "
          << sk_name << std::endl;}

        file.close();

      }

      if ( sk_name.size() > DFTBP_LC - 2) {std::cerr << "ERROR IN DFTB: SK FILENAME " << sk_name
        <<" IS TOO LONG " << std::endl;}

      for (int str_i = 0; str_i <  DFTBP_LC; str_i++){
        if (str_i < sk_name.size()) _dftb_options.skNames[str_i + counter * DFTBP_LC] = sk_name[str_i];
        else _dftb_options.skNames[str_i +  counter * DFTBP_LC] = DFTBP_PADCHAR;
      }
      counter++;
    }
  }
  _dftb_options.skNames[n_files * DFTBP_LC - 1] = '\0';

  // SPECIES NAMES
  _dftb_options.specieNames = (char *) malloc( _atomistic_structure->get_N_types() * DFTBP_MC * sizeof(char));

  counter = 0;
  for (int i = 0; i < _atomistic_structure->get_N_types(); i++){

    for (int str_i = 0; str_i <  DFTBP_MC - 1; str_i++){
      if (str_i <   _dftb_options.specieNameStrings[i].size() ) _dftb_options.specieNames[str_i + counter * DFTBP_MC] =  _dftb_options.specieNameStrings[i][str_i];
      else _dftb_options.specieNames[str_i +  counter * DFTBP_MC] = DFTBP_PADCHAR;
    }
    counter++;
  }

  _dftb_options.specieNames[ _atomistic_structure->get_N_types() * DFTBP_MC - 1] = '\0';


};



void
Dftb::build_structure_options()
{


  _dftb_options.nAtom = _atomistic_structure->get_N_atoms();
  _dftb_options.nType = _atomistic_structure->get_N_types();

  _dftb_options.coords = new double[_dftb_options.nAtom * 3];

  std::vector<Atom> basis;
  basis = _atomistic_structure->get_structure_atoms();

  double bohr_amstrong = Constants::bohr_radius * 1e10;
  //! Setting coordinates in DFTB format, coordinates are in atomic units
  for (int i = 0; i < _dftb_options.nAtom; i++){

    _dftb_options.coords[ (i*3) ] = basis[i].get_position(0) / bohr_amstrong;
    _dftb_options.coords[ (i*3) + 1 ] = basis[i].get_position(1) / bohr_amstrong;
    _dftb_options.coords[ (i*3) + 2 ] = basis[i].get_position(2) / bohr_amstrong;

  }

  //! Setting species in DFTB format
  _dftb_options.species = new int [ _dftb_options.nAtom ];

  for (int i = 0; i < _dftb_options.nAtom; i++){
    _dftb_options.species[i] = _atomistic_structure->get_type_index(basis[i].get_specie());

  }

  if (_atomistic_structure->is_periodic()) _dftb_options.iPeriodic = 1;
  else _dftb_options.iPeriodic = 0;

  //If system is not periodical lattice vectors are set to 0 by default in AtomisticStructure.
  //Check if huge value is needed instead of zero value
  _dftb_options.latVecs = _atomistic_structure->get_periodicity_vectors();
  //Convert latVecs in atomic units
  for (unsigned int i = 0; i < 9; i++) _dftb_options.latVecs[i] = _dftb_options.latVecs[i] / bohr_amstrong;

  //If lattice vectors are found in input file (in order x1, y1, z1, x2...), overwrite them
  if ( get_options().find_option("latvecs") ) {

    //If lattice vectors are specified, structure is considered periodical
    std::vector<double> tmp_latvecs;
    get_options().get_option("latvecs", tmp_latvecs);
    if (tmp_latvecs.size() != 9) std::cerr << "ERROR: latvecs must be a 9 doubles array" << std::endl;
    for (int i = 0; i < 9; i++) {_dftb_options.latVecs[i] = tmp_latvecs[i];}

  }


};


void
Dftb::build_input_options()
{
  std::cout << "build_input_options() begin " <<std::endl;
  std::string kpoints_path;
  std::ifstream file;
  std::string line, record;
  int counter = 0;
  double tmp;

  // Temperature is internally expressed in atomic units (Eh/Kb)
  // In input file it should be expressed in Kelvin
  tmp = get_options().get_option("eTemp", 300);
  _dftb_options.eTemp = tmp * (Constants::kb / Constants::Hartree );

  _dftb_options.orbResolved = get_options().get_option("orbResolved", 0);

  _dftb_options.skInterp = get_options().get_option("skInterp", 2);

  if (get_options().find_option("iPeriodic")){
    _dftb_options.iPeriodic = get_options().get_option("iPeriodic", 0);
    std::cout << "iPeriodic is " << _dftb_options.iPeriodic << std::endl;
  }

  //Maximum angolar momentum assigned in accordance with infos.data file
  _dftb_options.mAngs = new int[_dftb_options.nType];
  for (unsigned int i = 0; i < _dftb_options.nType; i++)   _dftb_options.mAngs[i] = _shell[_dftb_options.specieNameStrings[i]];


  if (_dftb_options.iPeriodic == 1){

    if (get_options().find_option("supercellfolding")) {

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

    else
      {

        //TODO: set a default k points set!!!!!!!!
        //KPOINTS are readen from an external file
        //Syntax is
        //<number of points>
        //value[1, x]  value[1,y]  value[1,z] weight
        // value[2,x]  value[2,y]  value[2,z] weight
        //.....

        kpoints_path = get_options().get_option("kpoints_path", "./kpoints.dat");

        file.open(kpoints_path.c_str());

        if ( !(file.is_open()) )
          {
            std::cerr << "COULD NOT FIND K POINTS FILE, PERFORMING GAMMA POINT"
            "CALCULATION " << kpoints_path << std::endl;
            file.close();

            _dftb_options.nkPoints = 1;
            _dftb_options.kPoints = new double [   _dftb_options.nkPoints * 3 ];
            _dftb_options.kWeights = new double[  _dftb_options.nkPoints];

            _dftb_options.kPoints[0] = 0.0;
            _dftb_options.kPoints[1] = 0.0;
            _dftb_options.kPoints[2] = 0.0;
            _dftb_options.kWeights[0] = 1.0;

          }

        else
          {

            getline(file, line);

            _dftb_options.nkPoints = atoi(line.c_str());

            if ( (_dftb_options.nkPoints <=0) )
              {
                std::cerr << "ERROR IN DFTB: NUMBER OF K POINTS IS ZERO OR NEGATIVE" << std::endl;
                exit(1);
              }

            _dftb_options.kPoints = new double [   _dftb_options.nkPoints * 3 ];
            _dftb_options.kWeights = new double[  _dftb_options.nkPoints];

            for (unsigned int i = 0; i < _dftb_options.nkPoints; i++)
              {
                getline(file,line);

                std::stringstream linestream(line);

                linestream >> _dftb_options.kPoints[0 +  counter * 3];

                linestream >> _dftb_options.kPoints[1 +  counter * 3];

                linestream >> _dftb_options.kPoints[2 +  counter * 3];

                linestream >> _dftb_options.kWeights[0 + counter];
                std::cout << "k point is " << _dftb_options.kWeights[0 + counter] <<
                _dftb_options.kPoints[0 +  counter * 3] << _dftb_options.kPoints[1 +  counter * 3] <<  _dftb_options.kPoints[2 +  counter * 3];
                std::cout << " and k points are " << _dftb_options.nkPoints << std::endl;
                if ( (_dftb_options.nkPoints <=0) )
                  {
                    std::cerr << "ERROR IN DFTB: NUMBER OF K POINTS IS ZERO OR NEGATIVE" << std::endl;
                    exit(1);
                  }

                counter = counter + 1;

              }

            //TODO: implementa KLines similmente a quanto avviene in DFTB+
          }

      }

  }
  std::cout << "build_input_options() done " <<std::endl;
};


void
Dftb::read_kpoints(void)
{

}


void
Dftb::obtain_hubbard_parameters(void)
{
  std::map<Shell, double> tmp_map;
  double* tmp_u_hub = NULL;

  tmp_u_hub = new double[_max_shell * _dftb_options.nType];
  inst->gethubbards(_dftb_options.nType, _max_shell, tmp_u_hub);

  //Note that order is due to the fact that in Fortran 2D arrays are stored by column
  for (unsigned int j = 0; j < _dftb_options.nType; j++)
    {
      tmp_map.clear();
      for (int i = 0; i < _max_shell; i++)
        {
          std::cout <<  _dftb_options.specieNameStrings[j] << std::endl;


          std::pair<std::string, std::map<Shell, double> > tmp_pair2;
          Shell tmp_shell;

          //Change integer in Shell type. Implicit conversion from int to enum is not allowed
          if ( i == 0 ) tmp_shell = S;
          else if ( i == 1 ) tmp_shell = P;
          else if ( i == 2) tmp_shell = D;
          else tmp_shell = NONE;

          tmp_map[tmp_shell] = * (tmp_u_hub + i + j);

        }
      _u_hub[Specie::string_to_specie[_dftb_options.specieNameStrings[j]]] = tmp_map;
    }

  delete[] tmp_u_hub; tmp_u_hub = NULL;

}


void
Dftb::print_dftb_options(void)
{
  std::cout << "DFTB_OPTIONS: " << std::endl;

  std::cout << "coords are " << std::endl;
  for (int i = 0; i < _dftb_options.nAtom * 3; i++){
    std::cout << "coords["<<i<<"] is " <<  _dftb_options.coords[ i ] << std::endl;
  }

  std::cout << "species are " << std::endl;
  for (int i = 0; i < _dftb_options.nAtom; i++){
    std::cout << "species["<<i<<"] is " <<  _dftb_options.species[ i ] << std::endl;
  }
  std::cout << "mAngs are " << std::endl;
  for (int i = 0; i < _dftb_options.nType; i++){
    std::cout << "mAngs["<<i<<"] is " <<  _dftb_options.mAngs[ i ] << std::endl;
  }

  int n_files = 0;

  std::cout << "latVecs are " << std::endl;
  for (int i = 0; i < 9; i++){
    std::cout << "latvecs[" <<i<<"]"<< _dftb_options.latVecs[i] << std::endl;
  }

  std::cout << "iPeriodic is " << _dftb_options.iPeriodic << std::endl;
  std::cout << "eTemp is " << _dftb_options.eTemp << std::endl;
  std::cout << "nAtom is " << _dftb_options.nAtom << std::endl;
  std::cout << "nType is " << _dftb_options.nType << std::endl;
  std::cout << "skInterp is " << _dftb_options.skInterp << std::endl;
  std::cout << "_max_shell is " << _max_shell << std::endl;

};


void
Dftb::add_shifts(void)
{
  double* shift_pnt = NULL;

  project_potential(_dftbp_coupling_options.potential_sim, "point");

  shift_pnt = new double[_pot_shift.size()];

  for (unsigned int i = 0; i < _pot_shift.size(); i++)
    {
      shift_pnt[i] = _pot_shift[i];
    }

  inst->setexternalshift(_dftb_options.nAtom, shift_pnt);

  delete[] shift_pnt;

}


#endif // ENABLE_DFTB

