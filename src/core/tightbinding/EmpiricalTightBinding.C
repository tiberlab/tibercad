
#include "tiber_config.h"


#ifdef ENABLE_UPTIGHT

#include "EmpiricalTightBinding.h"
#include "BoundaryProperties.h"
#include "PhysicalModel.h"
#include "SimulationOptions.h"
#include "EtbModel.h"
#include "SimulationOptions.h"
#include "UptWrapper.h"
#include "uptight.h"
#include "mesh.h"

#include <fstream>
#include <sstream>
#include <utility>
using namespace std;


//--------------------------------------------------------------

ETB::ETB(void)
: _upt_options()
{
  inst = new UptWrapper;
};


ETB::~ETB(void)
{
  delete inst;
  inst = NULL;
};

ETB::UptOptions::UptOptions(void)
:verbose(10),
 max_TB_order(2),
 harrison_flag(1),
 relat_flag(0),
 potential_flag(0),
 opt_flag(0),
 poldir(1)
{
  database_path = new char[UPT_LC];
  work_path = new char[UPT_LC];
  upt_filename = new char[UPT_MC];
  gen_outfile = new char[UPT_MC];
  c_axis = new double[3];
} 

ETB::UptOptions::~UptOptions(void)
{
  delete[] c_axis;
  delete[] database_path;
  delete[] work_path;
  delete[] upt_filename;
  delete[] gen_outfile;
};

ETB::UptSolverOptions::UptSolverOptions(void)
  :solver("upt lanczos"),
   n_vb(0),
   n_cb(0),
   min_iter(2),
   long_iter(30),
   max_iter(100000),
   guess_vb(0.0),
   guess_cb(0.0),
   fast_tol(1e-7),
   long_tol(1e-10),
   ort_tol(1e-4)
{
};

ETB::UptSolverOptions::~UptSolverOptions(void)
{
   
};

//-------------------------------------------------------------------------
PhysicalModel*
ETB::create_physical_model (const ModelOptions &options,
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

  std::cerr << "Empirical TB Initialisation" << std::endl;

  //_mesh = & ( get_environment().get_device().get_mesh());

  // Getting reference to atomistic structure for calculation
  get_atomistic_structure();

  // Setting options for DFTB+ calls
  parse_options();

  //  Initialize Dftb instance it and sets parameters
  inst->fill_param(_upt_options.verbose, _upt_options.database_path,
		   _upt_options.work_path, _upt_options.upt_filename, 
		   _upt_options.gen_outfile, _upt_options.max_TB_order,
		   _upt_options.harrison_flag, _upt_options.relat_flag, 
		   _upt_options.potential_flag, _upt_options.opt_flag,
		   _upt_options.poldir, _upt_options.c_axis);

  std::cout << "fill parameter done" << std::endl;

  //inst->addskdata(_dftb_options.skNames, _dftb_options.mAngs,
  //    _dftb_options.orbResolved, _dftb_options.skInterp, _dftb_options.nType);

  //std::cout << "addskdata done" << std::endl;

  std::cout << "init uptight begins" << std::endl;

  inst->inituptight();

  std::cout << "init uptight done" << std::endl;

};

//-------------------------------------------------------------------------
void ETB::do_solve(void){

#ifdef DEBUG
  std::cout << "Calling ETB->do_solve() " << std::endl;
#endif


  //TODO: non so come funziona il settaggio delle grandezze da calcolare quando
  //queste siano richieste da altri e non siano prettamente dati di output

  if (_upt_solver_options.solver.compare("upt lanczos") == 0) {

    inst->lanczos_diag(_upt_solver_options.n_vb, _upt_solver_options.n_cb, 
		 _upt_solver_options.guess_vb, _upt_solver_options.guess_cb,
		 _upt_solver_options.min_iter, _upt_solver_options.long_iter, 
		 _upt_solver_options.max_iter, _upt_solver_options.fast_tol, 
		 _upt_solver_options.long_tol, _upt_solver_options.ort_tol);

  }

#ifdef DEBUG
  std::cout << "ETB->do_solve() done" << std::endl;
#endif

};

//-------------------------------------------------------------------------

void ETB::parse_options(void)
{
  
  build_input_options();
  

#ifdef DEBUG
  print_upt_options();
#endif
  
};

//-------------------------------------------------------------------------
void
ETB::build_input_options()
{
  
  std::cout << "build_input_options() begin " <<std::endl;

  _upt_options.verbose = get_options().get_option("verbose", 10);
  _upt_options.max_TB_order = get_options().get_option("max_TB_order", 2);  
  _upt_options.harrison_flag = get_options().get_option("Harrison_scaling", true);
  _upt_options.relat_flag = get_options().get_option("relativistic", false); 
  _upt_options.opt_flag = get_options().get_option("optical transitions", false);   
  _upt_options.poldir = get_options().get_option("", false);   


  _upt_solver_options.n_vb =  get_solver_options().get_option("num_valence_eigenvalues", 0);
  _upt_solver_options.n_cb =  get_solver_options().get_option("num_conduction_eigenvalues", 0);
  _upt_solver_options.min_iter =  get_solver_options().get_option("min_iter", 2);
  _upt_solver_options.long_iter =  get_solver_options().get_option("long_iter", 30); 
  _upt_solver_options.max_iter =  get_solver_options().get_option("max_iter", 100000);
  _upt_solver_options.guess_vb =  get_solver_options().get_option("guess_valence", 0.0);
  _upt_solver_options.guess_cb =  get_solver_options().get_option("guess_conduction", 0.0);
  _upt_solver_options.fast_tol =  get_solver_options().get_option("fast_tolerance", 1e-7);
  _upt_solver_options.long_tol =  get_solver_options().get_option("long_tolerance", 1e-10);
  _upt_solver_options.ort_tol =  get_solver_options().get_option("orthogonality_tolerance", 1e-4);

  std::cout << "build_input_options() done " <<std::endl;

};

//-------------------------------------------------------------------------

void
ETB::print_upt_options(void)
{

  std::cout << "UPTIGHT_OPTIONS: " << std::endl;

  int n_files = 0;

  std::cout << "verbose is " << _upt_options.verbose << std::endl;
  std::cout << "max TB order is " << _upt_options.max_TB_order << std::endl;
  std::cout << "harrison scaling is " << _upt_options.harrison_flag << std::endl;
  std::cout << "relativistic is " << _upt_options.relat_flag << std::endl;
  std::cout << "external potential is " << _upt_options.potential_flag << std::endl;
  std::cout << "optical transitions is " << _upt_options.opt_flag << std::endl;
  std::cout << "polarization is along " << _upt_options.poldir << std::endl;
  std::cout << "database path is " << _upt_options.database_path << std::endl;
  std::cout << "work path is " << _upt_options.work_path << std::endl;
  std::cout << "upt filename is " << _upt_options.upt_filename << std::endl;
  std::cout << "gen output is " << _upt_options.gen_outfile << std::endl;
  std::cout << "c-axis is " << _upt_options.c_axis[1] << " " 
                 	    << _upt_options.c_axis[2] << " "
	                    << _upt_options.c_axis[3] << std::endl; 


   std::cout << "n valence " << _upt_solver_options.n_vb << std::endl;
   std::cout << "n conduction " <<  _upt_solver_options.n_cb << std::endl;
    std::cout << "min inter " << _upt_solver_options.min_iter << std::endl;
   std::cout << "long iter " <<  _upt_solver_options.long_iter << std::endl;
   std::cout << "max iter " <<  _upt_solver_options.max_iter << std::endl;
   std::cout << "guess valence " <<  _upt_solver_options.guess_vb << std::endl;
   std::cout << "guess conduction " <<  _upt_solver_options.guess_cb << std::endl;
   std::cout << "fast tol " <<  _upt_solver_options.fast_tol << std::endl;
   std::cout << "long tol " <<  _upt_solver_options.long_tol << std::endl;
   std::cout << "orth tol " <<  _upt_solver_options.ort_tol << std::endl;

};

//-------------------------------------------------------------------------
void
ETB::read_kpoints(void)
{

}

#endif
