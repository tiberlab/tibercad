
#include "tiber_config.h"


#ifdef ENABLE_UPTIGHT

#include "EmpiricalTightBinding.h"
//#include "BoundaryProperties.h"
#include "PhysicalModel.h"
#include "EtbModel.h"
#include "SimulationOptions.h"
#include "UptWrapper.h"
#include "uptight.h"
#include "mesh.h"

#include <fstream>
#include <sstream>
#include <utility>
#include <map>
#include <algorithm>
#include "Material.h"
#include "Alloy.h"
//#include <complex>
using namespace std;


//--------------------------------------------------------------

ETB::ETB(void)
: _upt_options()
{
  inst = new UptWrapper;
}


ETB::~ETB(void)
{
  delete inst;
  inst = NULL;
  _el_atomic_charges.clear();
  _hl_atomic_charges.clear();
  _band_shift.clear();
  _map_ID_Evb.clear();
  _map_ID_Ecb.clear();
  _ion_num_orbitals.clear();
}

ETB::UptOptions::UptOptions(void)
:verbose(10),
 max_TB_order(2),
 harrison_flag(1),
 relat_flag(0),
 potential_flag(0),
 opt_flag(0),
 poldir(1)
{
  c_axis = new double[3];
  c_axis[0]=0.0; c_axis[1]=0.0; c_axis[2]=1.0;
  database_path = new char[UPT_LC]; memset(database_path, ' ', UPT_LC);
  work_path = new char[UPT_LC];     memset(work_path, ' ', UPT_LC);
  out_path = new char[UPT_LC];      memset(out_path, ' ', UPT_LC);
  upt_filename = new char[UPT_MC];  memset(upt_filename, ' ', UPT_MC);
  gen_outfile = new char[UPT_MC];   memset(gen_outfile, ' ', UPT_MC);
  sparse_fmt = new char[UPT_MC];    memset(sparse_fmt, ' ', UPT_MC);

}

ETB::UptOptions::~UptOptions(void)
{
  delete[] c_axis;
  delete[] database_path;
  delete[] work_path;
  delete[] upt_filename;
  delete[] gen_outfile;
  delete[] sparse_fmt;
}

ETB::UptSolverOptions::UptSolverOptions(void)
  :solver("upt_lanczos"),
   n_vb(0),
   n_cb(0),
   min_iter(2),
   long_iter(30),
   max_iter(100000),
   guess_vb(0.0),
   guess_cb(0.0),
   fast_tol(1e-1),
   long_tol(1e-10),
   ort_tol(1e-4)
{
}

ETB::UptSolverOptions::~UptSolverOptions(void)
{

}

//-------------------------------------------------------------------------
PhysicalModel*
ETB::create_physical_model(const ModelOptions &options,
			     const Material* mat) const throw (ModelErrorException)
{

      ETBModel* model = dynamic_cast<ETBModel*> ( PhysicalModelInterface::create("etb",options) );

      if (model == NULL)
        throw ModelErrorException("TightBinding: ETB physical model is not created" );

      return model;

}


//-------------------------------------------------------------------------
void
ETB::do_init(void){

  std::cerr << "(TC) Empirical TB Initialisation..." << std::endl;

  TightBinding::do_init(); //gets mesh and _atomistic_structure

  if(_atomistic_structure==NULL)
    throw InitFailedException("ETB: atomistic structure not created");

  get_band_edges();  //reads band edges from database

  parse_options();   

#ifdef DEBUG
  print_upt_options();
#endif


  // Get database path from Database class
  std::string database_path = Database::get_default_search_path();
  std::string work_path = ".";
  std::string gen_outfile = "out.gen";
  std::string out_path = get_control().get_output_dir();

  if (database_path.size() > UPT_LC)
           throw InitFailedException("ETB: database search path too long");


  database_path.copy(_upt_options.database_path, database_path.size() );
  work_path.copy(_upt_options.work_path, work_path.size() );
  gen_outfile.copy(_upt_options.gen_outfile, gen_outfile.size() );
  out_path.copy(_upt_options.out_path, out_path.size() );



  //inst->addskdata(_dftb_options.skNames, _dftb_options.mAngs,
  //    _dftb_options.orbResolved, _dftb_options.skInterp, _dftb_options.nType);

  //std::cout << "addskdata done" << std::endl;

  //std::cout << "(TC) init uptight begins" << std::endl;

  _init = 1;
  _assemble = 1;

}

//-------------------------------------------------------------------------
void ETB::reinit(void){


  std::cout << "(TC) clean uptight data container" << std::endl;

  inst->cleanuptight();

  // checks that the strain simulation, if specified has been done
  if(_upt_options.strain_sim != "no_sim")
  {
    if( ! get_control().find_simulation(_upt_options.strain_sim)->is_solved() )
    throw InitFailedException("Strain model has not been solved");
  }

  std::string upt_filename;
  // Getting reference to atomistic structure for calculation
  if (_upg_filename.compare("none") != 0)
  {
    upt_filename = _upg_filename;
  }
  else
  {

    upt_filename = _atomistic_structure->get_name() + ".upg";

    std::cerr << "printing structure " << upt_filename << std::endl;

    _atomistic_structure->print_upg(upt_filename, _upt_options.etb_dataset);

    std::cout << "Number of atoms: " <<_atomistic_structure->get_N_atoms() << std::endl;

    std::cout << "Numb. without H: " <<_atomistic_structure->get_N_without_H()
	      << std::endl;
  }

  upt_filename.copy(_upt_options.upt_filename, upt_filename.size() );

  // temporary hack to consider c-axis orientations
  // By default wurtzites c-axis is along z. If 1-d then it is along x:
  if (get_environment().get_device().get_mesh().mesh_dimension() == 1)
    { _upt_options.c_axis[0]= 1.0;
      _upt_options.c_axis[1]= 0.0;
      _upt_options.c_axis[2]= 0.0;
    }


  std::cout << "(TC) init uptight begins" << std::endl;

  //  Set parameters for Uptight instance
  inst->fill_param(_upt_options.verbose, _upt_options.database_path,
		   _upt_options.work_path, _upt_options.out_path,
		   _upt_options.upt_filename,
		   _upt_options.gen_outfile, _upt_options.sparse_fmt,
		   _upt_options.max_TB_order, _upt_options.harrison_flag,
		   _upt_options.relat_flag, _upt_options.potential_flag,
		   _upt_options.opt_flag, _upt_options.poldir,
		   _upt_options.c_axis, _upt_options.check_bondmap);

  std::cout << "(TC) fill parameter done" << std::endl;

  inst->inituptight();

  std::cout << "(TC) init uptight done" << std::endl;

  _ion_num_orbitals.resize(_atomistic_structure->get_N_atoms(), 0);


  std::cout << "(TC) set orbitals per atom" << std::endl;

  inst->get_ion_numorbitals(_ion_num_orbitals);


  _N_without_H = _atomistic_structure->get_N_without_H();

}


//-------------------------------------------------------------------------
void ETB::do_solve(void){

  ModelOptions options;

  std::cout << "Tight-Binding calculations" << std::endl;

  std::set<ID> IDs = _atomistic_structure->get_IDset();
  std::set<ID>::iterator reg;

  for (reg = IDs.begin(); reg != IDs.end(); reg++)
  {
    std::cerr<< "Valence= "<<_map_ID_Evb[*reg] << std::endl;
    std::cerr<< "Conduct= "<<_map_ID_Ecb[*reg] << std::endl;   
  }

  std::cerr << "Vb Max= " << _vb_shift << std::endl;  

  print_upt_options();

  if (_init) reinit();

  _init = 0;

  if (_assemble) assemble(options);

  _assemble = 0;

  if (_upt_solver_options.solver.compare("upt_lanczos") == 0) {

    std::cout << "(TC) solving using lanczos" << std::endl;

    inst->lanczos_diag(_upt_solver_options.n_vb, _upt_solver_options.n_cb,
		 _upt_solver_options.guess_vb, _upt_solver_options.guess_cb,
		 _upt_solver_options.min_iter, _upt_solver_options.long_iter,
		 _upt_solver_options.max_iter, _upt_solver_options.fast_tol,
		 _upt_solver_options.long_tol, _upt_solver_options.ort_tol);

  }

  if (_upt_solver_options.solver.compare("read_old") == 0) {

    std::cout << "(TC) reading old states" << std::endl;

    inst->set_num_states(_upt_solver_options.n_vb, _upt_solver_options.n_cb);
    inst->read_old_states();

    int num_ev = _upt_solver_options.n_vb + _upt_solver_options.n_cb;
    std::complex<double> matel;

    for(int i=1; i<= num_ev; i++)
    {
	matel = inst->get_matel(i,i);
	std::cout << "eigval " << i << "= " << matel << std::endl;
    }
  }

#ifdef DEBUG
  std::cout << "(TC) ETB->do_solve() done" << std::endl;
  std::cout << "(TC) Copy solutions into _solutions" << std::endl;
#endif

  std::cerr<<"Pot min= "<<_pot_min<<std::endl;

  int hdim = inst->get_H_dim();
  int num_vb = _upt_solver_options.n_vb;
  int num_ev = _upt_solver_options.n_vb + _upt_solver_options.n_cb;
  double *eigvals = new double[num_ev];
  double *eigvects_re = new double[hdim*num_ev];
  double *eigvects_im = new double[hdim*num_ev];
  double *eigtmp_re = eigvects_re;
  double *eigtmp_im = eigvects_im;

  inst->get_states(num_ev,hdim,eigvals,eigvects_re,eigvects_im);

  _solution.resize(num_ev);

  for(int i=0; i< _upt_solver_options.n_vb; i++)
  {
    _solution[i].particle = "hl";
    _solution[i].statistics = "Fermi";
    _solution[i].eigen_energy = eigvals[i] + _vb_shift - _pot_min;
    _solution[i].eigen_vector.resize(hdim);
    _solution[i].temperature = _upt_options.temperature;

    if (i!= 0)
    {
      eigtmp_re += hdim;
      eigtmp_im += hdim;
    }

    for(int j=0; j<hdim;j++)
    {
      _solution[i].eigen_vector[j] = complex<double>(*(eigtmp_re+j),*(eigtmp_im+j));
    }

    if( (_upt_options.potential_flag) && !(get_options().find_option("hl_qfermi_level")))
    {
      _solution[i].electro_chem_pot = calculate_fermi_averaged(i);
    }
    else
    {
      _solution[i].electro_chem_pot = _upt_options.hl_chem_pot;
    }



  }


  for(int i=num_vb; i< num_vb + _upt_solver_options.n_cb; i++)
  {
    _solution[i].particle = "el";
    _solution[i].statistics = "Fermi";
    _solution[i].eigen_energy = eigvals[i] + _vb_shift - _pot_min;
    _solution[i].eigen_vector.resize(hdim);
    _solution[i].temperature = _upt_options.temperature;

    eigtmp_re += hdim;
    eigtmp_im += hdim;


    for(int j=0; j<hdim;j++)
    {
      _solution[i].eigen_vector[j] = complex<double>(*(eigtmp_re+j),*(eigtmp_im+j));

    }

    if( (_upt_options.potential_flag) && !(get_options().find_option("el_qfermi_level")))
    {
      _solution[i].electro_chem_pot = calculate_fermi_averaged(i);
    }
    else
    {
      _solution[i].electro_chem_pot =  _upt_options.el_chem_pot;
    }

  }

  // write state infos on screen.
  write_states();

  //Calculate electron and holes charge density on atoms (hydrogen not included)
  std::cerr << "Calculating quantum charge " << std::endl;

  _el_atomic_charges.resize(_N_without_H, 0.0);
  _hl_atomic_charges.resize(_N_without_H, 0.0);

  compute_atomic_charges("el", _el_atomic_charges);
  compute_atomic_charges("hl", _hl_atomic_charges);
  std::cerr << "done " << std::endl;

  //Print for debug charges on atoms
  double* charges;
  charges = new double[_atomistic_structure->get_N_atoms()];
  for (unsigned int i = 0; i < _N_without_H; i++) charges[i] = _el_atomic_charges[i] * 1000;
  for (unsigned int i = _N_without_H + 1; i < _atomistic_structure->get_N_atoms(); i++) charges[i] = 0.0;

  _atomistic_structure->print_structure("charges_el.xyz", charges);
  for (unsigned int i = 0; i < _N_without_H; i++) charges[i] = _hl_atomic_charges[i] * 1000;
  for (unsigned int i = _N_without_H + 1; i < _atomistic_structure->get_N_atoms(); i++) charges[i] = 0.0;
  _atomistic_structure->print_structure("charges_hl.xyz", charges);

  delete eigvals;
  delete eigvects_re;
  delete eigvects_im;

}

//-------------------------------------------------------------------------
void ETB::assemble(const ModelOptions& options)
{

  if( options.get_option("P_matrix",false) )
  {
    int poldir = options.get_option("poldir",0);
    inst->set_verbose(0);
    inst->compute_P_matrix(poldir);
    inst->set_verbose(_upt_options.verbose);
  }
  else
  {
    inst->clear_potential();

    if(_upt_options.band_shift_flag)
    {
      std::cerr<< "ETB: including band shifts" <<std::endl;
      add_band_shifts();
    }

    if(_upt_options.potential_flag)
    {
      std::cerr<< "ETB: passing potential" <<std::endl;
      add_pot_shifts();
    }
    inst->compute_H();

  }

}
//-------------------------------------------------------------------------
std::complex<double> ETB::calculate_matrix_element(const std::string& i_particle,
						   unsigned int i,
						   const std::string& j_particle,
						   unsigned int j)
{
  i++; j++; //this is done to set base vector to 1 (fortran mode)

  // the following operation is done to offset the electron states.
  // in the ETB solver-lib the hole states are stored first, then the electron states
  if (i_particle == "el" || i_particle == "electron")
  {
    i = i + _upt_solver_options.n_vb;
  }
  if (j_particle == "el" || j_particle == "electron")
  {
    j = j + _upt_solver_options.n_vb;
  }

  //std::cerr << "(" << i << "-" << j <<")"<<std::flush;

  return inst->get_matel(i,j);
}
//-------------------------------------------------------------------------
void ETB::solve_for_particle(const std::string& particle)
{

  if (particle == "el" || particle == "electron")
  {
    int temp_n_vb = _upt_solver_options.n_vb;
    _upt_solver_options.n_vb = 0;

    this->do_solve();

    _upt_solver_options.n_vb = temp_n_vb;

  }
  if (particle == "hl" || particle == "hole")
  {
    int temp_n_cb = _upt_solver_options.n_cb;
    _upt_solver_options.n_cb = 0;

    this->do_solve();

    _upt_solver_options.n_cb = temp_n_cb;
  }



}

//-------------------------------------------------------------------------
void ETB::do_plot(void){

  std::cout << "tightbinding do_plot()" << std::endl;
  TightBinding::do_plot();
  std::cout << "tightbinding do_plot()" << std::endl;

#ifdef DEBUG
  std::cout << "(TC) Calling ETB->do_plot() " << std::endl;
#endif

  const std::set<string>& plots = get_control().get_plotvariables();

  std::string outdir =get_environment().get_device()
                                       .get_control().get_output_dir();

  std::string file_name = outdir + "/states.dat";
  //write eigenvalues and population infos

  this->write_states(file_name);

  // write states in upg and cub formats
  inst->write_states();

}

//-------------------------------------------------------------------------

void ETB::parse_options(void)
{
  std::cout << "(TC) parse_options() begin...";

  _upg_filename = get_options().get_option("upg_filename", "none");

  _upt_options.verbose = get_options().get_option("verbose", verbose());

  _upt_options.etb_dataset = get_options().get_option("dataset","");
  _upt_options.max_TB_order = get_options().get_option("max_TB_order", 2);
  std::string sparse_fmt = get_options().get_option("sparse_format", "upper");
  sparse_fmt.copy(_upt_options.sparse_fmt, sparse_fmt.size() );

  _upt_options.check_bondmap = get_options().get_option("check_bondmap", false);
  _upt_options.harrison_flag = get_options().get_option("Harrison_scaling", true);
  _upt_options.relat_flag = get_options().get_option("relativistic", false);

  _upt_options.temperature = get_options().get_option("temperature",
						      SimulationOptions::temperature );

  //_upt_options.opt_flag = get_options().get_option("optical_transitions", false);
  //_upt_options.poldir = get_options().get_option("polarization_direction", 1);
  _upt_options.opt_flag = false; // these are set via OpticsTB
  _upt_options.poldir = 1;       //   "    "   "   "    "

  if (get_options().find_option("potential_simulation"))
  {
    _upt_options.potential_sim = get_options().get_option("potential_simulation","no_sim");
    _upt_options.potential_flag = true;
  }

  _upt_options.hl_chem_pot = get_options().get_option("hl_qfermi_level", 0.0);
  _upt_options.el_chem_pot = get_options().get_option("el_qfermi_level", 0.0);

  _upt_options.strain_sim = get_options().get_option("strain_model_name", "no_sim");

  // Solver options: "upt_lanczos", "read_old"
  _upt_solver_options.solver = get_solver_options().get_option("solver", "upt_lanczos");

  _upt_solver_options.n_vb =  get_solver_options().get_option("num_valence_eigenvalues", 0);
  if( _upt_solver_options.n_vb == 0) {
    _upt_solver_options.n_vb =  get_solver_options().get_option("num_hole_states", 0);
  }
  _upt_solver_options.n_cb =  get_solver_options().get_option("num_conduction_eigenvalues", 0);
  if( _upt_solver_options.n_cb == 0) {
    _upt_solver_options.n_cb =  get_solver_options().get_option("num_electron_states", 0);
  }
  _upt_solver_options.min_iter =  get_solver_options().get_option("min_iter", 2);
  _upt_solver_options.long_iter =  get_solver_options().get_option("long_iter", 30);
  _upt_solver_options.max_iter =  get_solver_options().get_option("max_iter", 100000);

  //---------------------------------------------------------------------------------------
  //computes educated guesses for valence and conduction bands edges
  std::set<ID> IDs = _atomistic_structure->get_IDset();
  std::set<ID>::iterator reg;

  double vb_max = -1000.0;
  for (reg = IDs.begin(); reg != IDs.end(); reg++)
  {
    if(_map_ID_Evb[*reg] > vb_max) vb_max = _map_ID_Evb[*reg];
  }
  double cb_min = 1000.0;
  for (reg = IDs.begin(); reg != IDs.end(); reg++)
  {
    if(_map_ID_Ecb[*reg] < cb_min) cb_min = _map_ID_Ecb[*reg];
  }
  
  // now vb_shift corresponds to maximum valence band edge
  _vb_shift = vb_max;
  
  if(cb_min<vb_max)
  {
    std::cerr<<"WARNING: bands overlap; cannot find good guess"<<std::endl;
  }
  else
  {
    vb_max = 1.0*(cb_min - vb_max)/5.0;
    cb_min = 4.0*vb_max;    
  }
  //---------------------------------------------------------------------------------------

  _upt_options.band_shift_flag = get_options().get_option("add_band_shifts", true);
  _upt_solver_options.guess_vb = get_solver_options().get_option("guess_valence", vb_max);
  _upt_solver_options.guess_cb = get_solver_options().get_option("guess_conduction", cb_min);

  // da togliere e leggere dal database: shift della banda di valenza (che e` 0)
  //_upt_options.vb_shift = get_options().get_option("vb_shift", 0.0);

  _upt_solver_options.fast_tol =  get_solver_options().get_option("fast_tolerance", 1e-1);
  _upt_solver_options.long_tol =  get_solver_options().get_option("long_tolerance", 1e-10);
  _upt_solver_options.ort_tol =  get_solver_options().get_option("orthogonality_tolerance", 1e-5);

  //Get projection_lenght for quantum charge projection (nm)
  _upt_options.projection_lenght = ("projection_lenght", 5.0);

  std::cout << "done" << std::endl;

}

//-------------------------------------------------------------------------

void
ETB::print_upt_options(void)
{

  std::cout << "(TC) UPTIGHT_OPTIONS: " << std::endl;

  int n_files = 0;

  std::cout << "verbose: " << _upt_options.verbose << std::endl;
  std::cout << "max TB order: " << _upt_options.max_TB_order << std::endl;
  std::cout << "harrison scaling: " << _upt_options.harrison_flag << std::endl;
  std::cout << "relativistic: " << _upt_options.relat_flag << std::endl;
  std::cout << "external potential: " << _upt_options.potential_flag << std::endl;
  //std::cout << "optical transitions: " << _upt_options.opt_flag << std::endl;
  //std::cout << "polarization is along: " << _upt_options.poldir << std::endl;
  //std::cout << "database path is " << _upt_options.database_path << std::endl;
  //std::cout << "work path is " << _upt_options.work_path << std::endl;
  //std::cout << "upt filename is " << _upt_options.upt_filename << std::endl;
  //std::cout << "gen output is " << _upt_options.gen_outfile << std::endl;
  std::cout << "c-axis: " << _upt_options.c_axis[0] << " "
                 	    << _upt_options.c_axis[1] << " "
	                    << _upt_options.c_axis[2] << std::endl;


   std::cout << "n valence: " << _upt_solver_options.n_vb << std::endl;
   std::cout << "n conduction: " <<  _upt_solver_options.n_cb << std::endl;
    std::cout << "min inter: " << _upt_solver_options.min_iter << std::endl;
   std::cout << "long iter: " <<  _upt_solver_options.long_iter << std::endl;
   std::cout << "max iter: " <<  _upt_solver_options.max_iter << std::endl;
   std::cout << "guess valence: " <<  _upt_solver_options.guess_vb << std::endl;
   std::cout << "guess conduction: " <<  _upt_solver_options.guess_cb << std::endl;
   std::cout << "fast tol: " <<  _upt_solver_options.fast_tol << std::endl;
   std::cout << "long tol: " <<  _upt_solver_options.long_tol << std::endl;
   std::cout << "orth tol: " <<  _upt_solver_options.ort_tol << std::endl;

}

//-------------------------------------------------------------------------

void
ETB::add_pot_shifts(void)
{
  project_potential(_upt_options.potential_sim, "point");

  inst->add_potential(_pot_shift);
}
//-------------------------------------------------------------------------
void
ETB::add_band_shifts(void)
{
  _band_shift.clear();
  unsigned int N = _atomistic_structure->get_N_atoms();
  _band_shift.resize(N, 0.0);

  const std::vector< Atom >& atom = _atomistic_structure->get_structure_atoms();

  for (unsigned int i = 0; i < N; i++)
  {
    // the sign is inverted because in upt the potential is subtracted
    _band_shift[i]= - _map_ID_Evb[atom[i].get_region_ID()] + _vb_shift;

    //std::cerr << _band_shift[i] << std::endl;
  }

  inst->add_potential(_band_shift);
}

//-------------------------------------------------------------------------

double
ETB::calculate_fermi_averaged(unsigned int i)
{

  double av, sum, atom_sum;
  unsigned int k, j, k_at, N_atoms_wo_H;

  sum = 0.0; k = 0; k_at = 0;

  N_atoms_wo_H = _atomistic_structure->get_N_without_H();

  if(_solution[i].particle == "el" || _solution[i].particle == "electron")
  {

    for (j = 0; j < N_atoms_wo_H; j++)
    {
      atom_sum = 0.0;
      for (k = k_at; k < k_at + _ion_num_orbitals[j]; k++)
      {
	atom_sum += std::norm(_solution[i].eigen_vector[k]);
      }

      k_at = k;
      sum +=  _el_chem_pot[j]*atom_sum;
    }

  }

  if(_solution[i].particle == "hl" || _solution[i].particle == "hole")
  {

    for (j = 0; j < N_atoms_wo_H; j++)
    {
      atom_sum = 0.0;
      for (k = k_at; k < k_at + _ion_num_orbitals[j]; k++)
      {
	atom_sum += std::norm(_solution[i].eigen_vector[k]);
      }

      k_at = k;
      sum +=  _hl_chem_pot[j]*atom_sum;

    }

  }

  // debug cross check on vector size
  //std::cerr<<"(TC-debug) N. orbitals = "<< k_at <<endl;
  //std::cerr<<"(TC-debug) Matrix dim. = "<<inst->get_H_dim()<<endl;

  return sum;

}

void
ETB::read_kpoints(void)
{

}

void
ETB::compute_atomic_charges(const std::string& particle, std::vector<double>& qmat)
{
  double atom_sum;
  unsigned int i, k, j, k_at, N_atoms_wo_H;
  unsigned int n = _solution.size();

  k = 0; k_at = 0;

  N_atoms_wo_H = _atomistic_structure->get_N_without_H();

  if(qmat.size() < N_atoms_wo_H)
    throw InitFailedException("(INT. ERROR) array mismatch in compute_atomic_charges");

  for(j=0; j<N_atoms_wo_H; j++){qmat[j] = 0.0; }
  for(i = 0; i < n; i++)
  {
    if(_solution[i].particle == particle)
    {
      k = 0; k_at = 0;
      for (j = 0; j < N_atoms_wo_H; j++)
      {
	atom_sum = 0.0;
	for (k = k_at; k < k_at + _ion_num_orbitals[j]; k++)
	{
	    atom_sum += std::norm(_solution[i].eigen_vector[k]);
	}

	k_at = k;

	double pop = Fermi(_solution[i].eigen_energy, _solution[i].electro_chem_pot,
		       _solution[i].temperature);

	if ( (_solution[i].particle == "hl") || (_solution[i].particle == "hole") )
	  {
	    pop = 1 - pop;
	  }

	qmat[j] += pop * atom_sum;
      }

    }
  }

}

void ETB::get_band_edges(void)
{
  
  AtomisticStructure* as = _atomistic_structure;
  std::set<ID> IDs = as->get_IDset();

  for(std::set<ID>::iterator reg = IDs.begin(); reg != IDs.end(); reg++)
  {
      Material* mat = as->get_device()->get_material( (*reg) );

      if (mat->is_alloy())
      {
	Alloy* alloy = static_cast<const Alloy*>(mat);

	Material* matA = alloy->get_component_A();
        Material* matB = alloy->get_component_B();
	double x = alloy->get_molar_fraction();

	std::cerr << "Fraction: " << x<< std::endl; 

	Database& dbA = matA->get_database();
	dbA.set_section("valenceband");
	double vbA = dbA.get("E_v",0.0);
	dbA.set_section("bandgap");
	double Eg = min( dbA.get("Eg_G", 1e6), dbA.get("Eg_X", 1e6));
	Eg = min( dbA.get("Eg_L", 1e6), Eg);
	double cbA = vbA + Eg;

	std::cerr << "Eg_a: " << Eg << std::endl;

	Database& dbB = matB->get_database();
	dbB.set_section("valenceband");
	double vbB = dbB.get("E_v",0.0);
	dbB.set_section("bandgap");
	Eg = min( dbB.get("Eg_G", 1e6), dbB.get("Eg_X", 1e6));
	Eg = min( dbB.get("Eg_L", 1e6), Eg);	
	double cbB = vbB + Eg;

	std::cerr << "Eg_b: " << Eg << std::endl;

	std::cerr << "Vb_a: " << vbA << std::endl; 
	std::cerr << "Vb_b: " << vbB << std::endl; 
	std::cerr << "Cb_a: " << cbA << std::endl; 
	std::cerr << "Cb_b: " << cbB << std::endl; 

	_map_ID_Evb[*reg] = x*vbA + (1-x)*vbB;
	_map_ID_Ecb[*reg] = x*cbA + (1-x)*cbB;
	
	dbA.set_section(""); dbB.set_section("");
      }
      else
      {
	Database& db = mat->get_database();
	db.set_section("valenceband");
	double vb = db.get("E_v",0.0);

	_map_ID_Evb[*reg] = vb;

	db.set_section("bandgap");
	double Eg = min( db.get("Eg_G", 1e6), db.get("Eg_X", 1e6));
	Eg = min( db.get("Eg_L", 1e6), Eg);

	std::cerr << "Eg: " << Eg << std::endl;
	std::cerr << "Vb: " << vb << std::endl; 
	std::cerr << "Cb: " << vb+Eg << std::endl;

	_map_ID_Ecb[*reg] = vb + Eg; // Gap at 0 K.
	
	db.set_section("");
      }

  }
  
}

ID
ETB::convert_variable_name_to_id(const std:: string& variable_name) const
{

  ID id = INVALID_ID;


  if (variable_name == "elDensity" )
    id  = EL_CH;
  if (variable_name == "hlDensity")
    id = HL_CH;


  return id;
}


//void
//ETB::get_solution_secure(const Elem* elem,
//    const std::set<ID>& ids, std::vector<std::map<ID, double> >& values)
//{
//
//  Point p = elem->centroid();
//  std::vector<Point> points(1);
//  points[0] = p;
//
//  get_solution_secure(elem,points,ids,values);
//
//}
//
//
//void
//ETB::get_solution_secure(const Elem* elem, const std::vector<Point>& p,
//    const std::set<ID>& ids, std::vector<std::map<ID, double> >& values)
//{
//  unsigned int np = p.size();
//  values.resize(np);
//
//  if (ids.count(EL_CH))
//    {
//      for (unsigned int n = 0; n < p.size(); n++)
//	{
//	  values[n][EL_CH] = build_rho("el",p[n]);
//	}
//    }
//
//  if (ids.count(HL_CH))
//    {
//      for (unsigned int n = 0; n < p.size(); n++)
//	{
//	  values[n][HL_CH] = build_rho("hl",p[n]);
//	}
//    }
//}


void
ETB::get_solution_secure(const Elem* elem,
    const std::set<ID>& ids, std::vector<std::map<ID, double> >& values)
{

  Point p = elem->centroid();
  std::vector<Point> points(1);
  points[0] = p;

  get_solution_secure(elem,points,ids,values);

}


void
ETB::get_solution_secure(const Elem* elem, const std::vector<Point>& p,
    const std::set<ID>& ids, std::vector<std::map<ID, double> >& values)
{
  unsigned int np = p.size();
  values.resize(np);


  if (ids.count(EL_CH))
    {
      for (unsigned int n = 0; n < p.size(); n++)
        {
          if (_elemental_result_el.find(elem) != _elemental_result_el.end())
            {
              values[n][EL_CH] = _elemental_result_el[elem];
            }
          else
            {
              _elemental_result_el[elem] = build_rho("el", elem->centroid());
              values[n][EL_CH] = _elemental_result_el[elem];
            }
        }
    }

  if (ids.count(HL_CH))
    {
      for (unsigned int n = 0; n < p.size(); n++)
        {
          if (_elemental_result_hl.find(elem) != _elemental_result_hl.end())
            {
              values[n][HL_CH] = _elemental_result_hl[elem];
            }
          else
            {
              _elemental_result_hl[elem] = build_rho("hl", elem->centroid());
              values[n][HL_CH] = _elemental_result_hl[elem];
            }

        }
    }
}


void
ETB::build_elemental_results(const std::set<std::string>& variables,
    std::vector<double>& results, std::vector<std::string>& legend)
{

std::cerr << "building elemental results " << std::endl;
std::set<std::string>::iterator mit;
for (mit=variables.begin() ; mit != variables.end(); mit++)
  {
    std::cout << "variable is " << *mit << std::endl;
  }
  // we only do something if we are on processor 0
  // TODO parallelize
  if (libMesh::processor_id() != 0)
    return;


  // if there is no mesh we can return immediately
  if (_mesh == NULL)
    return;



  unsigned int n_vars = 0;
  const unsigned int nn  = _mesh->n_active_elem();
  int el_ch = -1;
  int hl_ch = -1;

  if (variables.count("ElQuantumDensity"))
    {
      legend.resize( n_vars + 1 );
      legend[n_vars] = "ElQuantumDensity";
      el_ch = n_vars;
      n_vars++;
    }

  if (variables.count("HlQuantumDensity"))
      {
        legend.resize( n_vars + 1 );
        legend[n_vars] = "HlQuantumDensity";
        hl_ch = n_vars;
        n_vars++;
      }

  results.resize(nn * n_vars,0.0);

  MeshBase::const_element_iterator it =  _mesh->active_local_elements_begin();
  const MeshBase::const_element_iterator end =     _mesh->active_local_elements_end();


  unsigned int elem_number = 0;
  for ( ; it != end; ++it)
    {

      const Elem* elem = *it;

      unsigned int id = n_vars * elem_number;

      std::vector<std::map<ID, double> > values;

      std::set<ID> ids;
      ids.insert(EL_CH);
      ids.insert(HL_CH);

      //get_solution_secure(elem, ids, values);
      get_solution(elem, ids, values);

      //double charge = values[0][CHARGE];


      if (el_ch != -1)
        {
          //get_solution_secure(elem, ids, values);
          results[id + el_ch] = values[0][EL_CH];
        }

      if (hl_ch != -1)
             {
               //get_solution_secure(elem, ids, values);
               results[id + hl_ch] = values[0][HL_CH];
             }

      elem_number++;
    } //over element

  results.resize(elem_number * n_vars);

}


double
ETB::build_rho(const std::string& particle, const Point& r)
{
  double tau = _upt_options.projection_lenght / _atomistic_structure->get_scale();
  const double deltar_max = tau * 10; //Maximum cutoff distance
  double deltar, uhatom;
  double rho = 0.0;
  double x1, y1, z1;
  double x ,y, z;
  std::vector<double>* charges = NULL;

  if (particle == "el") charges = &(_el_atomic_charges);
  if (particle == "hl") charges = &(_hl_atomic_charges);


  //std::cerr << " tau is " << tau << " " ;

  x = r(0); y = r(1); z = r(2);

  if ((*charges).size() == 0)
    {
      std::cerr << "ERROR IN TIGHTBINDING: trying to build charge density "
      "but no mulliken charges are available" << std::endl;
      exit(1);
    }

  for (unsigned int iatm = 0; iatm  < _N_without_H; iatm++)
    {

      //std::cout << "rho before loop is " << rho << std::endl;
      //Getting Hubbard parameter
      //Up to now densities are mapped on orbital S
      //uhatom = _u_hub[_atomistic_structure->get_structure_atoms()[iatm].get_specie()][S];


      //Convert atom position to mesh units
      x1 = _atomistic_structure->get_structure_atoms()[iatm].get_position(1) / _atomistic_structure->get_scale();
      y1 = _atomistic_structure->get_structure_atoms()[iatm].get_position(2) / _atomistic_structure->get_scale();
      z1 = _atomistic_structure->get_structure_atoms()[iatm].get_position(3) / _atomistic_structure->get_scale();

      if (get_environment().get_device().get_mesh().mesh_dimension() == 1)
        {
          y1 = 0.0; z1 = 0.0;
        }
      if (get_environment().get_device().get_mesh().mesh_dimension() == 2)
              {
               z1 = 0.0;
              }
      //delta_r is already in mesh units in this way
      deltar = sqrt( (x - x1) * (x - x1) + (y - y1) * (y - y1) + (z - z1) * (z - z1));

      //Normalize to linear or superficial periodicity
      double normalize = 1.0;
      double* period = _atomistic_structure->get_periodicity_vectors();
      if (get_environment().get_device().get_mesh().mesh_dimension() == 1)
        {
          normalize = period[4]*period[8] - period[7]*period[5]
                    / _atomistic_structure->get_scale() / _atomistic_structure->get_scale();
        }
      if (get_environment().get_device().get_mesh().mesh_dimension() == 2)
              {
                normalize = period[8] / _atomistic_structure->get_scale();
              }


      //Also hubbard parameters (and tau) must be scaled in mesh units
      // (uhatom is in (atomic units)^(-1))
//      double tau = ( ( uhatom * ( 16.0 / 5.0 ) ) / (Constants::bohr_radius ) ) * get_control().get_device().get_mesh_units();

      if (deltar > deltar_max) continue;
      else
        {
          //std::cout << "uhatom is " << uhatom <<std::endl;
          rho = rho + ((*charges)[iatm] * tau * tau * tau * exp(-1.0 * tau * deltar));

        }
    }

  rho = rho / (8.0 * 3.141592653589793);

  //scale rho to q/cm^3 (carrier density)
  rho =  rho / (( get_control().get_device().get_mesh_units() * 100.0 ) *
        ( get_control().get_device().get_mesh_units() * 100.0 ) *
        ( get_control().get_device().get_mesh_units() * 100.0 ));

  return rho;

}


#endif
