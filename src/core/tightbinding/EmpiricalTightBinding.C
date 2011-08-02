
#include "tiber_config.h"


#ifdef ENABLE_UPTIGHT

#include "EmpiricalTightBinding.h"
#include "PhysicalModel.h"
#include "EtbModel.h"
#include "AtomisticStructure.h"
#include "SimulationOptions.h"
#include "TiberCad.h"
#include "UptWrapper.h"
#include "uptight.h"
#include "Material.h"
#include "Alloy.h"
#include "Messages.h"
#include "mesh.h"
#include "Macrostrain.h"


#include <fstream>
#include <sstream>
#include <utility>
#include <cstring>
#include <map>
#include <algorithm>
#include <limits>
 

//#include <complex>
using namespace std;


//--------------------------------------------------------------

ETB::ETB(const ModelOptions& options)
: TightBinding(options),
  _upt_options()
{
  inst = new UptWrapper;
}


ETB::~ETB(void)
{
  std::cout << "destructing";
  delete inst;
  inst = NULL;
  for (std::map<unsigned int, std::vector<double>>::iterator 
      i = _eigenvector_mag.begin(); i !=  _eigenvector_mag.end(); i++)
  {
    (*i).second.clear();
  }
  _eigenvector_mag.clear();
  _el_atomic_charges.clear();
  _hl_atomic_charges.clear();
  _band_shift.clear();
  _map_ID_Evb.clear();
  _map_ID_Ecb.clear();
  _ion_num_orbitals.clear();
  std::cout << "destructed" << std::endl;
}

ETB::UptOptions::UptOptions(void)
:verbose(10),
 max_TB_order(2),
 harrison_flag(1),
 relat_flag(0),
 potential_flag(0),
 opt_flag(0),
 poldir(1),
 hybrid_passivation(false),
 dg_scale(0.10),
 dg_onsite(-200.0),
 grid_step(0.5)
{
  //c_axis.reserve(3);
  c_axis[0]=0.0; c_axis[1]=0.0; c_axis[2]=1.0;
  //k_point.reserve(3);
  k_point[0]=0.0; k_point[1]=0.0; k_point[2]=0.0;
  database_path = new char[UPT_LC]; memset(database_path, UPT_PADCHAR, UPT_LC);
  work_path = new char[UPT_LC];     memset(work_path, UPT_PADCHAR, UPT_LC);
  out_path = new char[UPT_LC];      memset(out_path, UPT_PADCHAR, UPT_LC);
  upt_filename = new char[UPT_MC];  memset(upt_filename, UPT_PADCHAR, UPT_MC);
  gen_outfile = new char[UPT_MC];   memset(gen_outfile, UPT_PADCHAR, UPT_MC);
  sparse_fmt = new char[UPT_MC];    memset(sparse_fmt, UPT_PADCHAR, UPT_MC);
  //out_format = new char[UPT_MC];    memset(out_format, UPT_PADCHAR, UPT_MC);
}

ETB::UptOptions::~UptOptions(void)
{
  delete[] work_path;
  delete[] out_path;
  delete[] upt_filename;
  delete[] gen_outfile;
  delete[] sparse_fmt;
  //delete[] out_format;
  delete[] database_path;
}

ETB::UptSolverOptions::UptSolverOptions(void)
  :solver("upt_lanczos"),
   n_vb(0),
   n_cb(0),
   min_iter(2),
   long_iter(30),
   max_iter(15000),
   guess_vb(0.0),
   guess_cb(0.0),
   fast_tol(1e-1),
   long_tol(1e-10),
   ort_tol(1e-4),
   twice_cb(0),
   twice_vb(0)
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

  std::cerr << "(ETB) Empirical TB Initialisation..." << std::endl;

  //sanity check of complex number passing (this should be checked elsewere perhaps)
  Complex zz;
  double re,im;

  inst -> complex_test(re,im,zz);
  if (zz != Complex(re,im))
  {
    cerr<< zz <<  re << ", " << im << endl;
    throw InitFailedException("ETB: complex-passing test failed");
  }


  TightBinding::do_init(); //gets mesh and get_atomistic_structure()

  if(get_atomistic_structure()==NULL)
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
  std::string out_path = get_output_directory();

  if (database_path.size() > UPT_LC - 1)
           throw InitFailedException("ETB: database search path too long");


  std::size_t length = 0;
  length = database_path.copy(_upt_options.database_path, database_path.size() );
  //_upt_options.database_path[length] = '\0';
  length = work_path.copy(_upt_options.work_path, work_path.size() );
  //_upt_options.work_path[length] = '\0';
  length = gen_outfile.copy(_upt_options.gen_outfile, gen_outfile.size() );
  //_upt_options.gen_outfile[length] = '\0';
  length = out_path.copy(_upt_options.out_path, out_path.size() );
  //_upt_options.out_path[length] = '\0';

  _dim = get_mesh().mesh_dimension();

  //inst->addskdata(_dftb_options.skNames, _dftb_options.mAngs,
  //    _dftb_options.orbResolved, _dftb_options.skInterp, _dftb_options.nType);

  //std::cout << "addskdata done" << std::endl;

  //std::cout << "(ETB) init uptight begins" << std::endl;

  _init = 1;      // initialization must be called
  _assemble = 1;  // matrix assemble must be done

}

//-------------------------------------------------------------------------
void ETB::reinit(void){


  std::cout << "(ETB) clean uptight data container" << std::endl;

  inst->cleanuptight();

  _vb_shift = 0.0;

  // checks that the strain simulation, if specified has been done
  if(_upt_options.strain_sim != "no_sim")
  {
    SimulationInterface* strsim = find_simulation(_upt_options.strain_sim);
    if( (strsim != NULL) && (! strsim->is_solved()) )
      throw InitFailedException("Strain model has not been solved");
  }

  // Get the minimum CB and maximum VB edges
  ModelOptions solopts = get_solver_options();
  if (!solopts.find_option("guess_conduction") ||
      !solopts.find_option("guess_valence"))
  {
    double cb_min, vb_max;
    get_band_extrema(cb_min, vb_max);
    ostringstream os;
    os << "CB min = " << cb_min << "  VB max = " << vb_max << endl;
    Messages::info(os.str());

    if (!solopts.find_option("guess_conduction"))
      _upt_solver_options.guess_cb = cb_min;

    if (!solopts.find_option("guess_valence"))
      _upt_solver_options.guess_vb = vb_max;
  }

  std::string upt_filename;
  // Getting reference to atomistic structure for calculation
  if (_upg_filename.compare("none") != 0)
  {
    upt_filename = _upg_filename;
  }
  else
  {

    upt_filename = get_atomistic_structure()->get_name() + ".upg";

    Messages::info("printing structure "+upt_filename);

    get_atomistic_structure()->print_upg(upt_filename, _upt_options.etb_dataset);

    std::cout << "Number of atoms: " <<get_atomistic_structure()->get_N_atoms() << std::endl;

    std::cout << "Numb. without H: " <<get_atomistic_structure()->get_N_without_H()
	      << std::endl;
  }

  std::size_t length = 0;
  length = upt_filename.copy(_upt_options.upt_filename, upt_filename.size() );


  // Use rotated crystal to obtain direction of c-axis (z-direction)
  get_c_axis();

  std::cout << "(ETB) c-axis: "<<_upt_options.c_axis[0]<<" "
            <<_upt_options.c_axis[1]<<" "<<_upt_options.c_axis[2]<<std::endl;  
  
  std::cout << "(ETB) fill parameter " << std::endl;

  //  Set parameters for Uptight instance
  inst->fill_param(_upt_options.verbose, _upt_options.database_path,
		   _upt_options.work_path, _upt_options.out_path,
		   _upt_options.upt_filename,
		   _upt_options.gen_outfile, _upt_options.sparse_fmt,
		   _upt_options.max_TB_order, _upt_options.harrison_flag,
		   _upt_options.relat_flag, _upt_options.potential_flag,
		   _upt_options.opt_flag, _upt_options.poldir,
		   _upt_options.c_axis, _upt_options.check_bondmap,
                   _upt_options.dg_scale, _upt_options.dg_onsite,
                   _upt_options.hybrid_passivation);
  
	  
  inst->set_output((int) _upt_options.out_format, _upt_options.grid_step);

  //std::cout << "(ETB) fill parameter done" << std::endl;
  //std::cout.flush();

  if(inst->inituptight() != 0){ 
    throw InitFailedException("internal handlers do not match"); }


  std::cout << "(ETB) init uptight done" << std::endl;

  _ion_num_orbitals.resize(get_atomistic_structure()->get_N_atoms(), 0);

  std::cout << "(ETB) set orbitals per atom" << std::endl;

  inst->get_ion_numorbitals(_ion_num_orbitals);

  _N_without_H = get_atomistic_structure()->get_N_without_H();

  if (has_new_k())
  {
    _upt_options.k_point[0] = _k_vector[0];
    _upt_options.k_point[1] = _k_vector[1];
    _upt_options.k_point[2] = _k_vector[2];
  
    inst->set_kpoint(_upt_options.k_point); 
    k_is_old();
  }

  if(_upt_options.verbose>3) print_upt_options();

  _init = 0;
  _assemble = 1; //must reassemble matrix


}
//-------------------------------------------------------------------------
void ETB::get_c_axis(void)
{

  std::set<ID> IDs = get_atomistic_structure()->get_IDset();
  std::map<ID,Tensor1> cc;
  std::set<ID>::iterator reg = IDs.begin();

  for(; reg != IDs.end(); reg++)
  {
    const Material* mat = get_environment().get_device().get_material( (*reg) );
    Tensor2Gen RotM = (mat->get_rotated_crystal()).RotMatrix;
    Tensor1 zz(0); zz(3) = 1.0;
    cc[*reg] = RotM * zz; 
  }
  

  //reg--;
  //for(; reg != IDs.begin(); reg--)
  //{
  //std::set<ID>::iterator regp1 = reg; regp1++;
  //if (cc[*regp1](1) != cc[*reg](1) || cc[*regp1](2) != cc[*reg](2) || cc[*regp1](3) != cc[*reg](3) )
  //  throw InitFailedException("c-axis have different orientations in atomistic structure");
  //}  

  _upt_options.c_axis[0]= cc[*IDs.begin()](1);
  _upt_options.c_axis[1]= cc[*IDs.begin()](2);
  _upt_options.c_axis[2]= cc[*IDs.begin()](3);


}

//-------------------------------------------------------------------------
void ETB::do_solve(void){

  // check unused tags in solver_options
  const ModelOptions& sol_opt = get_solver_options();
  sol_opt.find_option("simulation"); // remove simulation name
  sol_opt.check_unused(); 
  

  ModelOptions options;

  std::cout << "(ETB) Tight-Binding calculations" << std::endl;

  std::set<ID> IDs = get_atomistic_structure()->get_IDset();
  std::set<ID>::iterator reg;

  //for (reg = IDs.begin(); reg != IDs.end(); reg++)
  //{
  //  std::cerr<< "Valence= "<<_map_ID_Evb[*reg] << std::endl;
  //  std::cerr<< "Conduct= "<<_map_ID_Ecb[*reg] << std::endl;   
  //}

  //std::cerr << "Vb Max= " << _vb_shift << std::endl;

  reinit(); 



  if (_assemble) assemble(options);

  if (_upt_solver_options.solver.compare("upt_lanczos") == 0) {

    std::cout << "(ETB) solving using lanczos" << std::endl;

    if (_upt_solver_options.twice_cb)
      std::cout << "(ETB) twice cb: true " << std::endl;

    if (_upt_solver_options.twice_vb)    
      std::cout << "(ETB) twice vb: true " << std::endl;

    inst->lanczos_diag(1, 1, _upt_solver_options.n_vb, _upt_solver_options.n_cb,
                     _upt_solver_options.guess_vb, _upt_solver_options.guess_cb,
                     _upt_solver_options.min_iter, _upt_solver_options.long_iter,
                     _upt_solver_options.max_iter, _upt_solver_options.fast_tol,
                     _upt_solver_options.long_tol, _upt_solver_options.ort_tol,
                     _upt_solver_options.twice_vb, _upt_solver_options.twice_cb);


  }

  if (_upt_solver_options.solver.compare("feast") == 0) 
    {
      Messages::info("Solving Tight Binding with FEAST eigensolver");
      inst->feast(_upt_solver_options.e_min, _upt_solver_options.e_max, _upt_solver_options.m0);
    } 

  if (_upt_solver_options.solver.compare("read_old") == 0) {

    std::cout << "(ETB) reading old states" << std::endl;

    inst->set_num_states(_upt_solver_options.n_vb, _upt_solver_options.n_cb);
    inst->read_old_states();

    int num_ev = _upt_solver_options.n_vb + _upt_solver_options.n_cb;
    Complex matel;

    for(int i=1; i<= num_ev; i++)
    {
	matel = inst->get_matel(i,i);
	std::cout << "eigval " << i << "= " << matel << std::endl;
    }
  }

#ifdef DEBUG
  std::cout << "(ETB) ETB->do_solve() done" << std::endl;
  std::cout << "(ETB) Copy solutions into _solutions" << std::endl;
#endif
  
  //std::cerr<<"Pot min= "<<_pot_min<<std::endl;

  int hdim = inst->get_H_dim();
  int num_vb = _upt_solver_options.n_vb;
  int num_ev = _upt_solver_options.n_vb + _upt_solver_options.n_cb;
  double *eigvals = new double[num_ev];
  Complex *eigvects = new Complex[hdim*num_ev];
  int *particles = new int[num_ev];
  Complex *eigtmp = eigvects;  // walking pointer on eigenvectors array

 
  inst->get_states(num_ev,hdim,eigvals,eigvects,particles);

  _solution.resize(num_ev);

  for(int i=0; i< num_ev; i++)
  {

   if(particles[i]==1)
   {
      _solution[i].particle = "el";
   }
   else if(particles[i]==-1) 
   {
      _solution[i].particle = "hl";
   }
   else
   {
      throw SolveFailedException("ETB: unkown particle type");
   }

    _solution[i].statistics = "Fermi";
    _solution[i].eigen_energy = eigvals[i];
    _solution[i].temperature = _upt_options.temperature;
    _solution[i].eigen_vector.resize(hdim);

    eigtmp = eigvects + i*hdim;

    for(int j=0; j<hdim;j++)
    {
      _solution[i].eigen_vector[j] = *(eigtmp+j);
    }

     
    if (_upt_options.potential_flag)
    {
       if(particles[i]==1 && !has_option("el_qfermi_level"))
       {
	  _solution[i].electro_chem_pot = calculate_fermi_averaged(i);
       }
       else
       {
          _solution[i].electro_chem_pot = _upt_options.el_chem_pot; 	
       } 
  
       if( particles[i]==-1 && !has_option("hl_qfermi_level"))
       {
          _solution[i].electro_chem_pot = calculate_fermi_averaged(i);
       }
       else
       {
          _solution[i].electro_chem_pot = _upt_options.hl_chem_pot; 
       }
     }
  }

  delete eigvals;
  delete eigvects;
  delete particles;

  //Set _solution_size in base class TightBinding 
  _solution_size = _solution.size();

  // write state infos on screen.
  write_states();

  //Calculate electron and holes charge density on atoms (hydrogen not included)
  _el_atomic_charges.resize(_N_without_H, 0.0);
  _hl_atomic_charges.resize(_N_without_H, 0.0);

  compute_atomic_charges("el", _el_atomic_charges);
  compute_atomic_charges("hl", _hl_atomic_charges);

  //_eigenvector_mag.insert(make_pair(1, std::vector<double>(_N_without_H)));
  for (unsigned int i = 0; i < _solution_size; i++)
  { 
    _eigenvector_mag.insert(make_pair(i, std::vector<double>(_N_without_H)));
    compute_eigenvector_mag(i, _eigenvector_mag[i]); 
  }


  //Print for debug charges on atoms
  double* charges;
  charges = new double[get_atomistic_structure()->get_N_atoms()];
  for (unsigned int i = 0; i < _N_without_H; i++) charges[i] = _el_atomic_charges[i];
  for (unsigned int i = _N_without_H + 1; i < get_atomistic_structure()->get_N_atoms(); i++) charges[i] = 0.0;

  get_atomistic_structure()->print_structure("charges_el.xyz", charges);
  for (unsigned int i = 0; i < _N_without_H; i++) charges[i] = _hl_atomic_charges[i];
  for (unsigned int i = _N_without_H + 1; i < get_atomistic_structure()->get_N_atoms(); i++) charges[i] = 0.0;
  get_atomistic_structure()->print_structure("charges_hl.xyz", charges);

}

//-------------------------------------------------------------------------
void ETB::do_assemble(const ModelOptions& options)
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
      std::cout<< "(ETB) including band shifts" <<std::endl;
      add_band_shifts();
    }

    if(_upt_options.potential_flag)
    {
      std::cout<< "(ETB) passing potential" <<std::endl;
      add_pot_shifts();
    }

    write_shifts();
    
    if(_upt_options.d_states_correction)
    {
      std::cout<< "(ETB) passing strain" <<std::endl;
      project_atom_strain();
    }


    inst->compute_H();

  }

  _assemble = 0;

}
//-------------------------------------------------------------------------
Complex ETB::calculate_matrix_element(const std::string& i_particle,
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
void
ETB::plot_globaldata(void)
{

  string outdir = get_output_directory();

  string filename(outdir + "/" + get_output_filename() + ".dat");
  ofstream file;
  file.open(filename.c_str());

  if (file.good())
  {
    // header
    file << "# TB eigenstates (" << get_name() << ")\n";

    file << "# Index" << setw(9)<< "Particle" << setw(12) << "EigenEnergy" 
         << setw(15) << "Occupation"
         << setw(12) << "FermiLevel" << setw(12) << "Temperature" << "\n";

    for (unsigned int i = 0; i < _solution.size(); i++)
    {
        file << setw(7) << i << setw(8) << _solution[i].particle 
             << setw(11) << _solution[i].eigen_energy << " "
             << setw(14) << get_population(i) << " "
             << setw(11) << _solution[i].electro_chem_pot << " "
             << setw(11) << _solution[i].temperature << "\n";
    }
  }

  Messages::info("(ETB) write wave functions on files");
  inst->write_states();

}

//-------------------------------------------------------------------------

void ETB::parse_options(void)
{
  //std::cout << "(ETB) parse_options() begin...";

  const ModelOptions& solopts = get_solver_options();
  

  _upg_filename = get_option("upg_filename", "none");

  _upt_options.verbose = get_option("verbose", SimulationOptions::verbose());

  _upt_options.etb_dataset = get_option("dataset","");
  _upt_options.max_TB_order = get_option("max_TB_order", 2);
  std::string sparse_fmt = get_option("sparse_format", "upper");
  sparse_fmt.copy(_upt_options.sparse_fmt, sparse_fmt.size() );

  _upt_options.check_bondmap = get_option("check_bondmap", false);

  _upt_options.relat_flag = get_option("relativistic", true);

  _upt_options.temperature = get_option("temperature",
			SimulationOptions::temperature );

  //_upt_options.opt_flag = get_option("optical_transitions", false);
  //_upt_options.poldir = get_option("polarization_direction", 1);
  _upt_options.opt_flag = false; // these are set via OpticsTB
  _upt_options.poldir = 1;       //   "    "   "   "    "

  if (has_option("potential_simulation"))
  {
    _upt_options.potential_sim = get_option("potential_simulation","");
    _upt_options.potential_flag = true;
  }

  _upt_options.hl_chem_pot = get_option("hl_qfermi_level", 0.0);
  _upt_options.el_chem_pot = get_option("el_qfermi_level", 0.0);

  _upt_options.strain_sim = get_option("strain_model_name", "no_sim");
  _upt_options.strain_sim = get_option("strain_simulation", _upt_options.strain_sim);

  if (_upt_options.strain_sim == "no_sim")
  {
    _upt_options.harrison_flag = get_option("Harrison_scaling", false);
    _upt_options.d_states_correction = get_option("d_splitting",false); 
  }
  else
  {
    _upt_options.harrison_flag = get_option("Harrison_scaling", true);
    _upt_options.d_states_correction = get_option("d_splitting",true);     
  }


  // Dangling bond scaling
  _upt_options.dg_scale = get_option("dangling_bond_scaling",100);
  _upt_options.dg_onsite = get_option("dangling_bond_onsite",-200.0);    

  //Choose passivation model
  std::string passivation_model = get_option("passivation_model","hydrogen");
  if ( passivation_model == "hybrid" )
      {_upt_options.hybrid_passivation = true;}

  // Solver options: "upt_lanczos", "read_old"
  _upt_solver_options.solver = solopts.get_option("solver", "upt_lanczos");

  _upt_solver_options.n_vb =  solopts.get_option("num_valence_eigenvalues", 0);
  if( _upt_solver_options.n_vb == 0) {
    _upt_solver_options.n_vb =  solopts.get_option("num_hole_states", 0);
  }
  _upt_solver_options.n_cb =  solopts.get_option("num_conduction_eigenvalues", 0);
  if( _upt_solver_options.n_cb == 0) {
    _upt_solver_options.n_cb =  solopts.get_option("num_electron_states", 0);
  }
  _upt_solver_options.min_iter =  solopts.get_option("min_iter", 2);
  _upt_solver_options.long_iter =  solopts.get_option("long_iter", 30);
  _upt_solver_options.max_iter =  solopts.get_option("max_iter", 15000);

  bool flag = solopts.get_option("remove_folded_sols_conduction",false);
  if (flag) _upt_solver_options.twice_cb = 1;
  flag = solopts.get_option("remove_folded_sols_valence",false);
  if (flag) _upt_solver_options.twice_vb = 1;

  //Feast options
  _upt_solver_options.e_min =  solopts.get_option("Emin", 0.0);
  _upt_solver_options.e_max =  solopts.get_option("Emax", 3.0);
  _upt_solver_options.m0 =  solopts.get_option("subspace", 100);
  //---------------------------------------------------------------------------------------
  // output wavevetors format
  std::string out_fmt = get_option("jmol_output_format", "jvxl");

  // the next lines are ludicrous in order to solve a crazy problem 
  // in passing char* to fortran  
  if (out_fmt=="jvxl") _upt_options.out_format = JVXL;
  else if (out_fmt=="cube") _upt_options.out_format = CUBE;
  else throw InitFailedException("cannot recognize format "+out_fmt);
  
  _upt_options.grid_step = get_option("jmol_grid_step", 0.5);


  //---------------------------------------------------------------------------------------
  //computes educated guesses for valence and conduction bands edges
  std::set<ID> IDs = get_atomistic_structure()->get_IDset();
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
  //_vb_shift = vb_max;
  _vb_shift = 0.0;
  
  if(cb_min<vb_max)
  {
    std::cerr<<"WARNING: bands overlap; cannot find good guess"<<std::endl;
  }
  //else
  //{
  //  vb_max = 1.0*(cb_min - vb_max)/5.0;
  //  cb_min = 4.0*vb_max;
  //}
  //---------------------------------------------------------------------------------------

  _upt_options.band_shift_flag = get_option("add_band_shifts", true);
  _upt_solver_options.guess_vb = solopts.get_option("guess_valence", vb_max);
  _upt_solver_options.guess_cb = solopts.get_option("guess_conduction", cb_min);

  // da togliere e leggere dal database: shift della banda di valenza (che e` 0)
  //_upt_options.vb_shift = options.get_option("vb_shift", 0.0);

  _upt_solver_options.fast_tol =  solopts.get_option("fast_tolerance", 1e-1);
  _upt_solver_options.long_tol =  solopts.get_option("long_tolerance", 1e-10);
  _upt_solver_options.ort_tol =  solopts.get_option("orthogonality_tolerance", 1e-5);

  //Get projection_length for quantum charge projection (nm)
  _upt_options.projection_length = get_option("projection_length", 5.0);

  //std::cout << "Projection lenght set to " <<  _upt_options.projection_length << std::endl;
  //std::cout << "done" << std::endl;

  // get kpoint
  RealVectorValue k_vec;
  get_option("k_vector",k_vec);
  set_k_point(k_vec);

  // get kpoint as parameter
  get_parameter("k_x", _upt_options.k_point[0], _k_vector[0] );
  get_parameter("k_y", _upt_options.k_point[1], _k_vector[1] );
  get_parameter("k_z", _upt_options.k_point[2], _k_vector[2] );


}

//-------------------------------------------------------------------------

void
ETB::print_upt_options(void)
{

  std::cout << "(ETB) UPTIGHT_OPTIONS: " << std::endl;

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
  unsigned int N = get_atomistic_structure()->get_N_atoms();
  _band_shift.resize(N, 0.0);

  const std::vector< Atom >& atom = get_atomistic_structure()->get_structure_atoms();

  const Bondmap& b_map = get_atomistic_structure()->get_bond_map();

  for (unsigned int i = 0; i < N; i++)
  {
   
    if (atom[i].get_specie() == H)
    {
      // the sign is inverted because in upt the potential is subtracted
      _band_shift[i]= - _map_ID_Evb[atom[b_map[i][0]].get_region_ID()] + _vb_shift;
    }
    else
    {
      _band_shift[i]= - _map_ID_Evb[atom[i].get_region_ID()] + _vb_shift;
    } 
  }
  inst->add_potential(_band_shift);
}

//-------------------------------------------------------------------------
void
ETB::write_shifts(void)
{

  std::string outdir = TiberCad::get_output_dir();
  std::string file_name = outdir + "/pot_on_atoms.dat";
  std::ofstream file;

  file.open(file_name.c_str());

  if(_upt_options.band_shift_flag && _upt_options.potential_flag)
  {

     for(int i=0; i < _band_shift.size(); i++)
     {
        file << _band_shift[i] + _pot_shift[i] << "  "
             << _band_shift[i] << "  "
             << _pot_shift[i] << std::endl;
     }    
  }
 
  if(_upt_options.band_shift_flag && !_upt_options.potential_flag)
  {
     for(int i=0; i < _band_shift.size(); i++)
     {
         file << _band_shift[i] << std::endl;
     }    
   }

     
   if(_upt_options.potential_flag && !_upt_options.band_shift_flag)
   {
      for(int i=0; i < _pot_shift.size(); i++)
      {
          file << _pot_shift[i] << std::endl;
      }    
   }

  file.close();

}

//-------------------------------------------------------------------------------------------/
void
ETB::project_atom_strain(void)
{  
  std::vector<Atom>& structure = get_atomistic_structure()->get_structure_atoms();

  Macrostrain* strsim = dynamic_cast<Macrostrain*> (find_simulation(_upt_options.strain_sim));

  
  unsigned int Number_of_atoms =  get_atomistic_structure()->get_N_without_H();

  std::vector<double> exx(Number_of_atoms, 0.0);
  std::vector<double> eyy(Number_of_atoms, 0.0);  
  std::vector<double> ezz(Number_of_atoms, 0.0);

  unsigned int k = 0; 

  for (unsigned int i = 0; i < structure.size() ; i++)
  { 
    
    if (structure[i].get_specie() != H && structure[i].get_elem() != NULL)
    {
      Tensor2Sym epsilon = strsim->get_strain(structure[i].get_elem());
      k++;

      exx[k] = epsilon(1,1);
      eyy[k] = epsilon(2,2);
      ezz[k] = epsilon(3,3);
    }
  
  }
  
  inst->set_strain(exx,eyy,ezz);
  

}
 

double
ETB::calculate_fermi_averaged(unsigned int i)
{

  double av, sum, atom_sum;
  unsigned int k, j, k_at, N_atoms_wo_H;

  sum = 0.0; k = 0; k_at = 0;

  N_atoms_wo_H = get_atomistic_structure()->get_N_without_H();

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

  N_atoms_wo_H = get_atomistic_structure()->get_N_without_H();

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

void
ETB::compute_eigenvector_mag(unsigned int eigenstate, std::vector<double>& densatm)
{
  double atom_sum;
  unsigned int i, k, j, k_at, N_atoms_wo_H;
  unsigned int n = _solution.size();

  k = 0; k_at = 0;

  N_atoms_wo_H = get_atomistic_structure()->get_N_without_H();

  if(densatm.size() < N_atoms_wo_H)
    throw InitFailedException("(INT. ERROR) array mismatch in compute_state_density");

  if (eigenstate >= n)
    throw InitFailedException("Eigenstate index is larger than number of available eigenstates");

  for(j=0; j<N_atoms_wo_H; j++){densatm[j] = 0.0; }
      k = 0; k_at = 0;
      for (j = 0; j < N_atoms_wo_H; j++)
      {
	atom_sum = 0.0;
	for (k = k_at; k < k_at + _ion_num_orbitals[j]; k++)
	{
	    atom_sum += std::norm(_solution[eigenstate].eigen_vector[k]);
	}

	k_at = k;

	densatm[j] += atom_sum;
      }
}

void ETB::get_band_edges(void)
{
  
  AtomisticStructure* as = get_atomistic_structure();
  std::set<ID> IDs = as->get_IDset();

  for(std::set<ID>::iterator reg = IDs.begin(); reg != IDs.end(); reg++)
  {
      const Material* mat = as->get_device()->get_material( (*reg) );

      /*if (mat->is_alloy())
      {
	const Alloy* alloy = static_cast<const Alloy*>(mat);

	const Material* matA = alloy->get_component_A();
        const Material* matB = alloy->get_component_B();
	double x = alloy->get_molar_fraction();

	std::cerr << "Fraction: " << x<< std::endl; 

	Database dbA = matA->get_database();
	dbA.set_section("valenceband");
	double vbA = dbA.get("E_v",0.0);
	dbA.set_section("bandgap");
	double Eg = min( dbA.get("Eg_G", 1e6), dbA.get("Eg_X", 1e6));
	Eg = min( dbA.get("Eg_L", 1e6), Eg);
	double cbA = vbA + Eg;

	std::cerr << "Eg_a: " << Eg << std::endl;

	Database dbB = matB->get_database();
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
      else*/
      {
	Database db = mat->get_database();
	db.set_section("valenceband");
	double vb = db.get("E_v",0.0);

	_map_ID_Evb[*reg] = vb;

	db.set_section("bandgap");
	double Eg = min( db.get("Eg_G", 1e6), db.get("Eg_X", 1e6));
	Eg = min( db.get("Eg_L", 1e6), Eg);

	//std::cerr << "Eg: " << Eg << std::endl;
	//std::cerr << "Vb: " << vb << std::endl; 
	//std::cerr << "Cb: " << vb+Eg << std::endl;

	_map_ID_Ecb[*reg] = vb + Eg; // Gap at 0 K.
	
	db.set_section("");
      }

  }
  
}



void
ETB::get_band_extrema(double& cb_min, double& vb_max)
{
  if (_upt_options.potential_flag)
  {
    SimulationInterface* sim =
        SimulationInterface::find_simulation(_upt_options.potential_sim);

    Messages::info("(ETB) get band extrema from "+_upt_options.potential_sim);

    string cbedge("Ec");
    string vbedge("Ev");

    // here we assume that the simulation exists (it was checked before)
    ID cb_id = sim->get_solution_id(cbedge);
    ID vb_id = sim->get_solution_id(vbedge);

    if ((cb_id == INVALID_ID) || (vb_id == INVALID_ID))
      throw RuntimeException("Simulation \'" + sim->get_name() +
          "\' lacks band edge solution variables.");

    vb_max = -numeric_limits<double>::max();
    cb_min = numeric_limits<double>::max();

    //map<ID, vector<double> > bandedges;
    vector<double> cb_edges(8,0.0);
    vector<double> vb_edges(8,0.0);

    MeshBase::const_element_iterator it = get_mesh().active_elements_begin();
    const MeshBase::const_element_iterator end = get_mesh().active_elements_end();

    for ( ; it != end; ++it)
    { 
      const Elem* elem = *it;
      vector<Point> p(elem->n_nodes());

      for (size_t i = 0; i < elem->n_nodes(); ++i)
        p[i] = elem->point(i);

      sim->get_solution(elem, vb_id, vb_edges, p);
      sim->get_solution(elem, cb_id, cb_edges, p);

      //for (size_t i = 0; i < elem->n_nodes(); ++i)
      //  std::cout<<vb_edges[i]<<" ";
      //std::cout<<std::endl;
      
      for (size_t i = 0; i < elem->n_nodes(); ++i)
      {
        double vb = vb_edges[i];
        double cb = cb_edges[i];
        vb_max = (vb > vb_max) ? vb : vb_max;
        cb_min = (cb < cb_min) ? cb : cb_min;
      }

      //std::cout<<vb_max<<std::endl;
    }
  }
  else
  {
    //computes guess from database band edges
    std::set<ID> IDs = get_atomistic_structure()->get_IDset();
    std::set<ID>::iterator reg;

    vb_max = -1000.0;
    for (reg = IDs.begin(); reg != IDs.end(); reg++)
    {
      if(_map_ID_Evb[*reg] > vb_max) vb_max = _map_ID_Evb[*reg];
    }

    cb_min = 1000.0;
    for (reg = IDs.begin(); reg != IDs.end(); reg++)
    {
      if(_map_ID_Ecb[*reg] < cb_min) cb_min = _map_ID_Ecb[*reg];
    }

  }

  if(cb_min < vb_max)
  {
    Messages::warning("Bands overlap, cannot find good guess fo ETB");
  }
}




void
ETB::get_solution_secure(const Elem* elem,
    std::map<ID, std::vector<double> >& values,
    const std::vector<Point>& p)
{

  unsigned int np = p.size();

  //CELL values
  //
  if (values.count(ElQuantumDensity))
    {
	      if (_dim == 3)
                  values[ElQuantumDensity][0] = build_rho3d(_el_atomic_charges, elem->centroid());
              else if (_dim == 2)
                  values[ElQuantumDensity][0] = build_rho2d(_el_atomic_charges, elem->centroid());
              else if (_dim == 1)
                  values[ElQuantumDensity][0] = build_average_rho1d(_el_atomic_charges, elem);
              // probably unnecessary, otherwise should throw exception
              //else std::cerr << "Unknown number of dimensions, values not assigned in" 
              //  " ETB::get_solution_secure" << std::endl;
    }

  if (values.count(HlQuantumDensity))
    {
              if (_dim == 3)
                values[HlQuantumDensity][0] = build_rho3d(_hl_atomic_charges, elem->centroid());
              else if (_dim == 2)
                values[HlQuantumDensity][0] = build_rho2d(_hl_atomic_charges, elem->centroid());
              else if (_dim == 1)
                values[HlQuantumDensity][0] = build_average_rho1d(_hl_atomic_charges, elem);
              else std::cerr <<"Unknown number of dimensions, values not assigned in" 
                " ETG::get_solutionsecure" << std::endl;
    }
 
  if (values.count(MeshStates))
  {
    for (unsigned int i = 0; i < _solution.size(); i++)
    {
              if (_dim == 3)
                values[MeshStates][i] = build_rho3d(_eigenvector_mag[i], elem->centroid());
              else if (_dim == 2)
                values[MeshStates][i] = build_rho2d(_eigenvector_mag[i], elem->centroid());
              else if (_dim == 1)
                values[MeshStates][i] = build_average_rho1d(_eigenvector_mag[i], elem);
              else std::cerr <<"Unknown number of dimensions, values not assigned in" 
                " ETG::get_solutionsecure" << std::endl;
    }
  }

}


double
ETB::build_rho3d(const std::vector<double>& tb_density, const Point& r)
{
  double tau = 1.0 / (_upt_options.projection_length / get_atomistic_structure()->get_scale());
  const double deltar_max = tau * 10; //Maximum cutoff distance
  double deltar, uhatom;
  double rho = 0.0;
  double x1, y1, z1;
  double x ,y, z;
  
 
  x = r(0); y = r(1); z = r(2);

  if (tb_density.size() != _N_without_H)
    {
      std::cerr << "ERROR IN ETB: trying to build mesh density "
      "but tb density has wrong size" << std::endl;
      exit(1);
    }

  for (unsigned int iatm = 0; iatm  < _N_without_H; iatm++)
    {
      //Convert atom position to mesh units
      x1 = get_atomistic_structure()->get_structure_atoms()[iatm].get_position(1) / get_atomistic_structure()->get_scale();
      y1 = get_atomistic_structure()->get_structure_atoms()[iatm].get_position(2) / get_atomistic_structure()->get_scale();
      z1 = get_atomistic_structure()->get_structure_atoms()[iatm].get_position(3) / get_atomistic_structure()->get_scale();

      //delta_r is already in mesh units in this way
      deltar = sqrt( (x - x1) * (x - x1) + (y - y1) * (y - y1) + (z - z1) * (z - z1));

      //Also hubbard parameters (and tau) must be scaled in mesh units
      // (uhatom is in (atomic units)^(-1))

      if (deltar > deltar_max) continue;
      else
        {
          //std::cout << "uhatom is " << uhatom <<std::endl;
          rho = rho + (tb_density[iatm] * tau * tau * tau * exp(-1.0 * deltar * tau));

        }
    }

  rho = rho / (8.0 * 3.141592653589793);

  //scale rho to q/cm^3 (carrier density)
  double mesh_units = 100.0 * get_mesh_units();
  rho =  rho / ( mesh_units * mesh_units * mesh_units );
  return rho;

}


double
ETB::build_rho2d(const std::vector<double>& tb_density, const Point& r)
{
  double tau = 1.0 / (_upt_options.projection_length / get_atomistic_structure()->get_scale());
  const double deltar_max = tau * 10; //Maximum cutoff distance
  double deltar, uhatom;
  double rho = 0.0;
  double x1, y1, z1;
  double x ,y, z;


  x = r(0); y = r(1); z = 0.0;

  if (tb_density.size() != _N_without_H)
    {
      std::cerr << "ERROR IN ETB: trying to mesh density "
      "but tb_density has wrong size" << std::endl;
      exit(1);
    }

  for (unsigned int iatm = 0; iatm  < _N_without_H; iatm++)
    {
      //Convert atom position to mesh units
      x1 = get_atomistic_structure()->get_structure_atoms()[iatm].get_position(1) / get_atomistic_structure()->get_scale();
      y1 = get_atomistic_structure()->get_structure_atoms()[iatm].get_position(2) / get_atomistic_structure()->get_scale();
      z1 = 0.0;

      //delta_r is already in mesh units in this way
      deltar = sqrt( (x - x1) * (x - x1) + (y - y1) * (y - y1) );


      if (deltar > deltar_max) continue;
      else
        {
          //std::cout << "uhatom is " << uhatom <<std::endl;
          rho = rho + (tb_density[iatm] * tau * tau * exp(-1.0 * deltar * tau));

        }
    }

  rho = rho / (2.0 * 3.141592653589793);

  //scale rho to q/cm^3 (carrier density)
  double mesh_units = 100.0 * get_mesh_units();
  rho =  rho / ( mesh_units * mesh_units * mesh_units );

  return rho;

}


double
ETB::build_average_rho1d(const std::vector<double>& tb_density, const Elem* elem)
{
  //Projection length in mesh units
  double proj_length = (_upt_options.projection_length / get_atomistic_structure()->get_scale());
  double tau = 1.0 / proj_length;
  const double deltar_max = proj_length * 10; //Maximum cutoff distance
  double deltar, uhatom, l;
  double rho = 0.0;
  double x_atm, x1, x2;

  assert(elem->n_nodes() == 2);

  x1 = elem->point(0)(0); x2 = elem->point(1)(0);

  //just order nodes in ascending order, using "x_atm" temporary
  x_atm = x1;
  if (x1 > x2) {x1 = x2; x2 = x_atm;}
  l = x2 - x1;

  if (tb_density.size() != _N_without_H)
    {
      std::cerr << "ERROR IN ETB: trying to build mesh density "
      "but tb_density has wrong size" << std::endl;
      exit(1);
    }

  double* period = get_atomistic_structure()->get_periodicity_vectors();

  for (unsigned int iatm = 0; iatm  < _N_without_H; iatm++)
    {

      //Convert atom position to mesh units
      x_atm = get_atomistic_structure()->get_structure_atoms()[iatm].get_position(1) / get_atomistic_structure()->get_scale();

      //delta_r is already in mesh units in this way
      deltar = min( abs( (x_atm - x1) ), abs( (x_atm - x2) ) );

      if (deltar > deltar_max) continue;
      else
      {
        //Calculate charge
        //1D formulation is the solution of (DQ*tau/2)*int_{x1}^{x2}{exp(-tau*|x-x_atm|)dx}

        if (x_atm >= x2)
        {
          rho = rho + (0.5) * (tb_density[iatm]) * l * ( + exp(tau * (x2 - x_atm)) - exp(tau * (x1 - x_atm)) );
        }
        else if (x_atm <= x1)
        {
          rho = rho - (0.5) * (tb_density[iatm]) * l * ( + exp(-tau * (x2 - x_atm)) - exp(-tau * (x1 - x_atm)) );
        }
        else
        //x_atm between x1 and x2
        {
          rho = rho + (0.5) * (tb_density[iatm]) * l * (2.0 - exp(-tau * (x2 - x_atm)) - exp(tau * (x1 - x_atm)) );
        }
      }
    }

  double normalize = 1.0;
  //normalize = (period[4]*period[8] - period[7]*period[5]);

  //rho = rho / normalize;

  //double mesh_units = 100.0 * get_mesh_units();
  //rho =  rho / ( mesh_units * mesh_units * mesh_units );
  //std::cout << "I builded average 1d" << std::endl;
  return rho;

}


void
ETB::do_setup_solution_variables(void)
{
  unsigned int dim = get_mesh().mesh_dimension();
  string units("/cm");
  if (dim == 2)
     units = "/cm^2";
  else if (dim == 3)
     units = "/cm^3";

  declare_solution(ElQuantumDensity, REAL, CELL, "q"+units);
  declare_solution(HlQuantumDensity, REAL, CELL, "q"+units);
  declare_solution(MeshStates, NTUPLE, CELL, "1"+units, 1);

  if (plot_solution("ProbabilityDensity"))
	    add_plot_variable(MeshStates);

}

#endif

