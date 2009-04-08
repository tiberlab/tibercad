
#include "tiber_config.h"


#ifdef ENABLE_UPTIGHT

#include "EmpiricalTightBinding.h"
//#include "BoundaryProperties.h"
#include "PhysicalModel.h"
#include "EtbModel.h"
#include "SimulationOptions.h"
#include "UptWrapper.h"
#include "uptight.h"
//#include "mesh.h"

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

  std::cout << "(TC) parse options..."; 
  parse_options();
  std::cout << "done" << std::endl;

#ifdef DEBUG
  print_upt_options();
#endif  

  std::string upt_filename;
  // Getting reference to atomistic structure for calculation
  if (_upg_filename.compare("none") != 0)
  {
    upt_filename = _upg_filename;
  }
  else
  {
    get_atomistic_structure();
    upt_filename = _atomistic_structure->get_name() + ".upg";
  }

  // Get mesh informations
  //_mesh = & ( get_environment().get_device().get_mesh());

  // Get database path from Database class
  std::string database_path = Database::get_default_search_path();
  std::string work_path = ".";
  std::string gen_outfile = "out.gen";

  if (database_path.size() > UPT_LC) 
           throw InitFailedException("ETB: database search path too long");

  
  database_path.copy(_upt_options.database_path, database_path.size() ); 
  work_path.copy(_upt_options.work_path, work_path.size() );  
  gen_outfile.copy(_upt_options.gen_outfile, gen_outfile.size() ); 
  upt_filename.copy(_upt_options.upt_filename, upt_filename.size() );

  //  Set parameters for Uptight instance
  inst->fill_param(_upt_options.verbose, _upt_options.database_path, 
		   _upt_options.work_path, _upt_options.upt_filename, 
		   _upt_options.gen_outfile, _upt_options.sparse_fmt, 
		   _upt_options.max_TB_order, _upt_options.harrison_flag, 
		   _upt_options.relat_flag, _upt_options.potential_flag, 
		   _upt_options.opt_flag, _upt_options.poldir, 
		   _upt_options.c_axis, _upt_options.check_bondmap);

  std::cout << "(TC) fill parameter done" << std::endl;

  //inst->addskdata(_dftb_options.skNames, _dftb_options.mAngs,
  //    _dftb_options.orbResolved, _dftb_options.skInterp, _dftb_options.nType);

  //std::cout << "addskdata done" << std::endl;

  std::cout << "(TC) init uptight begins" << std::endl;

  inst->inituptight();

  std::cout << "(TC) init uptight done" << std::endl;

  _init = 0;
  _assemble = 1;
  
}

//-------------------------------------------------------------------------
void ETB::reinit(void){
  
  inst->cleanuptight(); 
  inst->inituptight();

}


//-------------------------------------------------------------------------
void ETB::do_solve(void){

  ModelOptions options;

  if (_assemble) assemble(options);

  _assemble = 0;
 

#ifdef DEBUG
  std::cout << "(TC) Calling ETB->do_solve() " << std::endl;
#endif

  //TODO: non so come funziona il settaggio delle grandezze da calcolare quando
  //queste siano richieste da altri e non siano prettamente dati di output

  if (_upt_solver_options.solver.compare("upt_lanczos") == 0) {

    inst->lanczos_diag(_upt_solver_options.n_vb, _upt_solver_options.n_cb, 
		 _upt_solver_options.guess_vb, _upt_solver_options.guess_cb,
		 _upt_solver_options.min_iter, _upt_solver_options.long_iter, 
		 _upt_solver_options.max_iter, _upt_solver_options.fast_tol, 
		 _upt_solver_options.long_tol, _upt_solver_options.ort_tol);

  }

#ifdef DEBUG
  std::cout << "(TC) ETB->do_solve() done" << std::endl;
  std::cout << "(TC) Copy solutions into _solutions" << std::endl;
#endif

  int hdim = inst->get_H_dim();
  int num_ev = _upt_solver_options.n_vb + _upt_solver_options.n_cb;
  double *eigvals = new double[num_ev];
  double *eigvects_re = new double[hdim*num_ev];
  double *eigvects_im = new double[hdim*num_ev];
  
  inst->get_states(num_ev,hdim,eigvals,eigvects_re,eigvects_im);

  //_solution.reserve(num_ev);

  for(int i=0; i< _upt_solver_options.n_vb; i++)
  {  
    _solution[i].particle = "hl";
    _solution[i].statistics = "Fermi";
    _solution[i].eigen_energy = eigvals[i];
    _solution[i].eigen_vector.reserve(hdim);
    for(int j=0; j<hdim;j++)
    {
      _solution[i].eigen_vector[j] = complex<double>(eigvects_re[j],eigvects_im[j]);  
    }
  }

  for(int i=0; i< _upt_solver_options.n_cb; i++)
  {  
    _solution[i].particle = "el";
    _solution[i].statistics = "Fermi";
    _solution[i].eigen_energy = eigvals[i];
    _solution[i].eigen_vector.reserve(hdim);
    for(int j=0; j<hdim;j++)
    {
      _solution[i].eigen_vector[j] = complex<double>(eigvects_re[j],eigvects_im[j]);
    
    }
  }
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
    inst->compute_P_matrix(poldir);  
  }
  else
  {
    if(_upt_options.potential_flag) add_shifts();

    inst->compute_H();
  }

}
//-------------------------------------------------------------------------
std::complex<double> ETB::calculate_matrix_element(const std::string& i_particle,
						   unsigned int i, 
						   const std::string& j_particle,
						   unsigned int j)
{
  double re,im;
  std::complex<double> matel;

  if (i_particle == "el" || i_particle == "electron")
  {
    i = i + _upt_solver_options.n_vb;
  }
  if (j_particle == "el" || j_particle == "electron")
  {
    j = j + _upt_solver_options.n_vb;
  }  

  inst->get_matel(i,j,re,im);
    
  matel = complex<double>(re,im);
  
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

  if (plots.find("tbstates") != plots.end())
  {
    inst->write_states();
  }
}

//-------------------------------------------------------------------------

void ETB::parse_options(void)
{
     std::cout << "(TC) parse_options() begin...";

  _upg_filename = get_options().get_option("upg_filename", "none"); 

  _upt_options.verbose = get_options().get_option("verbose", 10);
  _upt_options.max_TB_order = get_options().get_option("max_TB_order", 2); 
  std::string sparse_fmt = get_options().get_option("sparse_format", "upper");
  sparse_fmt.copy(_upt_options.sparse_fmt, sparse_fmt.size() ); 

  _upt_options.check_bondmap = get_options().get_option("check_bondmap", false);
  _upt_options.harrison_flag = get_options().get_option("Harrison_scaling", true);
  _upt_options.relat_flag = get_options().get_option("relativistic", false); 
  //_upt_options.opt_flag = get_options().get_option("optical_transitions", false);   
  //_upt_options.poldir = get_options().get_option("polarization_direction", 1);   
  _upt_options.opt_flag = false; // these are set via OpticsTB
  _upt_options.poldir = 1;       //   "    "   "   "    "

  if (get_options().get_option("potential_simulation", false))
  {
    _upt_options.potential_sim = get_options().get_option("potential_simulation","no_sim");
    _upt_options.potential_flag = true;
  }

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

  shift_pnt = new double[_pot_shift.size()];

  for (unsigned int i = 0; i < _pot_shift.size(); i++)
    {
      shift_pnt[i] = _pot_shift[i];
    }

  inst->add_potential(_pot_shift.size(), shift_pnt);

  delete[] shift_pnt;

}


void
ETB::read_kpoints(void)
{

}



#endif
