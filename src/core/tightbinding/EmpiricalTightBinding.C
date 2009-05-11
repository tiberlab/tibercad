
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

  // make sure an atomistic structure do exist
  get_atomistic_structure();

  if(_atomistic_structure==NULL)
    throw InitFailedException("ETB: atomistic structure not created");

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

    _atomistic_structure->print_structure(upt_filename);
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

  inst->get_ion_numorbitals(_ion_num_orbitals);

}


//-------------------------------------------------------------------------
void ETB::do_solve(void){

  ModelOptions options;

  std::cout << "Tight-Binding calculations" << std::endl;

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
    _solution[i].eigen_energy = eigvals[i] + _upt_options.vb_shift - _pot_min;
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

    if(_upt_options.potential_flag)
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
    _solution[i].eigen_energy = eigvals[i] + _upt_options.vb_shift - _pot_min;
    _solution[i].eigen_vector.resize(hdim);
    _solution[i].temperature = _upt_options.temperature;

    eigtmp_re += hdim;
    eigtmp_im += hdim;


    for(int j=0; j<hdim;j++)
    {
      _solution[i].eigen_vector[j] = complex<double>(*(eigtmp_re+j),*(eigtmp_im+j));

    }

    if(_upt_options.potential_flag)
    {
      _solution[i].electro_chem_pot = calculate_fermi_averaged(i);
    }
    else
    {
      _solution[i].electro_chem_pot =  _upt_options.el_chem_pot;
    }

  }

  write_states();

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
    if(_upt_options.potential_flag)
    {
      std::cerr<< "ETB: passing potential" <<std::endl;
      add_shifts();
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

#ifdef DEBUG
  std::cout << "(TC) Calling ETB->do_plot() " << std::endl;
#endif

  const std::set<string>& plots = get_control().get_plotvariables();

  inst->write_states();

}

//-------------------------------------------------------------------------

void ETB::parse_options(void)
{
  std::cout << "(TC) parse_options() begin...";

  _upg_filename = get_options().get_option("upg_filename", "none");

  _upt_options.verbose = get_options().get_option("verbose", 1);
  _upt_options.max_TB_order = get_options().get_option("max_TB_order", 2);
  std::string sparse_fmt = get_options().get_option("sparse_format", "upper");
  sparse_fmt.copy(_upt_options.sparse_fmt, sparse_fmt.size() );

  _upt_options.check_bondmap = get_options().get_option("check_bondmap", false);
  _upt_options.harrison_flag = get_options().get_option("Harrison_scaling", true);
  _upt_options.relat_flag = get_options().get_option("relativistic", false);
  // da togliere e leggere dal database: shift della banda di valenza (che e` 0)
  _upt_options.vb_shift = get_options().get_option("vb_shift", 0.0);

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
  _upt_solver_options.guess_vb =  get_solver_options().get_option("guess_valence", 0.0);
  _upt_solver_options.guess_cb =  get_solver_options().get_option("guess_conduction", 0.0);
  _upt_solver_options.fast_tol =  get_solver_options().get_option("fast_tolerance", 1e-1);
  _upt_solver_options.long_tol =  get_solver_options().get_option("long_tolerance", 1e-10);
  _upt_solver_options.ort_tol =  get_solver_options().get_option("orthogonality_tolerance", 1e-5);

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
  std::cout << "optical transitions: " << _upt_options.opt_flag << std::endl;
  std::cout << "polarization is along: " << _upt_options.poldir << std::endl;
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
ETB::add_shifts(void)
{
  double* shift_pnt = NULL;

  project_potential(_upt_options.potential_sim, "point");

  inst->add_potential(_pot_shift);

}

double
ETB::calculate_fermi_averaged(unsigned int i)
{

  double av, sum, atom_sum;
  unsigned int k, j, k_at;

  sum = 0.0; k = 0; k_at = 0;


  if(_solution[i].particle == "el" || _solution[i].particle == "electron")
  {

    for (j = 0; j < _el_chem_pot.size(); j++)
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

    for (j = 0; j < _hl_chem_pot.size(); j++)
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



#endif
