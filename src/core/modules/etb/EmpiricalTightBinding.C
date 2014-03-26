// $Id$

#include "EmpiricalTightBinding.h"
#include "PhysicalModel.h"
#include "AtomisticStructure.h"
#include "SimulationOptions.h"
#include "TiberCad.h"
#include "UptWrapper.h"
#include "uptight.h"
#include "Material.h"
#include "Alloy.h"
#include "Messages.h"
#include "mesh.h"
#include "EigenSolver.h"
#include "RotatedCrystal.h"

#include <fstream>
#include <sstream>
#include <utility>
#include <cstring>
#include <map>
#include <algorithm>
#include <limits>
 
#include "TiberModule.h"

//#include <complex>
using namespace std;


//--------------------------------------------------------------

ETB::ETB(const ModelOptions& options)
: TightBinding(options),
  _upt_options()
{
  inst = new UptWrapper;

  has_solution_vector(false);
}

ETB* ETB::create(const ModelOptions& options)
{
  return new ETB(options);
}

ETB::~ETB(void)
{
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
 dg_scale(1.0),
 dg_onsite(-200.0),
 grid_step(0.5)
{
  //c_axis.reserve(3);
  c_axis[0]=0.0; c_axis[1]=0.0; c_axis[2]=1.0;
  //k_point.reserve(3);
  k_point[0]=0.0; k_point[1]=0.0; k_point[2]=0.0;
  default_path = new char[UPT_LC];  memset(default_path, UPT_PADCHAR, UPT_LC);
  database_path = new char[UPT_LC]; memset(database_path, UPT_PADCHAR, UPT_LC);
  work_path = new char[UPT_LC];     memset(work_path, UPT_PADCHAR, UPT_LC);
  load_path = new char[UPT_LC];     memset(load_path, UPT_PADCHAR, UPT_LC);
  out_path = new char[UPT_LC];      memset(out_path, UPT_PADCHAR, UPT_LC);
  upt_filename = new char[UPT_MC];  memset(upt_filename, UPT_PADCHAR, UPT_MC);
  gen_outfile = new char[UPT_MC];   memset(gen_outfile, UPT_PADCHAR, UPT_MC);
  sparse_fmt = new char[UPT_MC];    memset(sparse_fmt, UPT_PADCHAR, UPT_MC);
  //out_format = new char[UPT_MC];    memset(out_format, UPT_PADCHAR, UPT_MC);
}

ETB::UptOptions::~UptOptions(void)
{
  delete[] work_path;
  delete[] load_path;
  delete[] out_path;
  delete[] upt_filename;
  delete[] gen_outfile;
  delete[] sparse_fmt;
  //delete[] out_format;
  delete[] database_path;
  delete[] default_path;
}

ETB::UptSolverOptions::UptSolverOptions(void)
  :solver("upt_lanczos"),
   start_vb(1),
   start_cb(1),
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
   twice_vb(0),
   dynamic(0)	
{
}

ETB::UptSolverOptions::~UptSolverOptions(void)
{

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


  TightBinding::do_init(); 

  if(get_atomistic_structure()==NULL)
    throw InitFailedException("ETB: atomistic structure not created");

  get_bulk_edges();  //reads band edges from database

  parse_options();   

#ifdef DEBUG
  print_upt_options();
#endif


  // Get database path from Database class
  //  std::string database_path = Database::get_default_search_path();
  std::string default_path = get_option("default_path",Database::get_default_search_path());
  std::string database_path = get_option("database_path",Database::get_search_path());
  std::string work_path = ".";
  std::string gen_outfile = "out.gen";
  std::string out_path = get_output_directory();

  if (database_path.size() > UPT_LC - 1)
           throw InitFailedException("ETB: database search path too long");


  std::size_t length = 0;
  length = default_path.copy(_upt_options.default_path, default_path.size() );
  length = database_path.copy(_upt_options.database_path, database_path.size() );
  length = work_path.copy(_upt_options.work_path, work_path.size() );
  length = gen_outfile.copy(_upt_options.gen_outfile, gen_outfile.size() );
  length = out_path.copy(_upt_options.out_path, out_path.size() );

  _dim = get_mesh().mesh_dimension();

  Messages::info("(ETB) creating map elem->atoms");

  build_map_elem_atoms(_upt_options.projection_length);

  

  _init = true;      // initialization must be called
  _assemble = true;  // matrix assemble must be done

}

//-------------------------------------------------------------------------
void ETB::do_reinit(void)
{

  std::cerr << "(ETB) clean uptight data container" << std::endl;

  inst->cleanuptight();

  if (!includes_regions(get_atomistic_structure()->get_IDset()))
  {
     std::cerr << "(ETB) restrict on active regions" << std::endl;

     std::set<ID> reg_ids;
     get_region_ids(reg_ids);

     get_atomistic_structure()->restrict(reg_ids);
    
     std::cerr<< "(ETB) Build Map Elem->Atom "<<std::endl;

     double Rmax = build_map_elem_atoms(_upt_options.projection_length);
     
     std::cerr<< "(ETB) Rmax= "<<Rmax<<std::endl;
  }

  _vb_shift = 0.0;

  // checks that the strain simulation, if specified has been done
  if(_upt_options.strain_sim != "no_sim")
  {
    SimulationInterface* strsim = _strain_int.get_simulation();

    if( ! strsim->is_solved() ) 
      throw InitFailedException("Strain model has not been solved");
  }

  // checks that the potential simulation, if specified has been done
  if(_upt_options.potential_sim != "no_sim")
  {
    if( ! _dd_int->is_solved() )
      throw InitFailedException(_upt_options.potential_sim+" model has not been solved");
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
      _upt_solver_options.guess_cb = cb_min - 0.3;

    if (!solopts.find_option("guess_valence"))
      _upt_solver_options.guess_vb = vb_max + 0.3;
  }

  std::string upt_filename;
  // Getting reference to atomistic structure for calculation
  if (_upg_filename.compare("none") != 0)
  {
    upt_filename = _upg_filename;
  }
  else
  {

    upt_filename = get_output_directory() + "/" +
        get_atomistic_structure()->get_name() + ".upg";

    //Messages::info("(ETB) printing structure "+upt_filename);

    get_atomistic_structure()->print_upg(upt_filename, _upt_options.etb_dataset, 
                                                      !_upt_options.band_shift_flag);

    std::cout << "(ETB) Number of atoms: " <<get_atomistic_structure()->get_N_atoms() << std::endl;

    std::cout << "(ETB) Numb. without H: " <<get_atomistic_structure()->get_N_without_H()
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
  inst->set_paths(_upt_options.default_path, _upt_options.database_path, 
                  _upt_options.work_path, _upt_options.out_path);
  
  inst->fill_param(_upt_options.verbose, _upt_options.upt_filename,
		   _upt_options.gen_outfile, _upt_options.sparse_fmt,
		   _upt_options.max_TB_order, _upt_options.harrison_flag,
		   _upt_options.relat_flag, _upt_options.potential_flag,
		   _upt_options.opt_flag, _upt_options.poldir,
		   _upt_options.c_axis, _upt_options.check_bondmap,
       _upt_options.dg_scale, _upt_options.dg_onsite,
       _upt_options.hybrid_passivation);
 
  std::cout << "(ETB) set solver flag "<< _upt_solver_options.solver_flag << std::endl;
  inst->set_solver_flag(_upt_solver_options.solver_flag); 
	  
  inst->set_output((int) _upt_options.out_format, _upt_options.grid_step);

  //std::cout << "(ETB) fill parameter done" << std::endl;
  //std::cout.flush();

  if(inst->inituptight() != 0){ 
    throw InitFailedException("internal handlers do not match"); }

  _ion_num_orbitals.resize(get_atomistic_structure()->get_N_atoms(), 0);

  inst->get_ion_numorbitals(_ion_num_orbitals);

  _N_without_H = get_atomistic_structure()->get_N_without_H();

  if (has_new_k())
  {
    Point kp(get_k_point(true));
    cerr << kp << endl;
    _upt_options.k_point[0] = kp(0);
    _upt_options.k_point[1] = kp(1);
    _upt_options.k_point[2] = kp(2);
  
    inst->set_kpoint(_upt_options.k_point); 
    // this is now private, and maybe not really needed here:
    //k_is_old();
  }

  if(_upt_options.verbose>3) print_upt_options();

  _init = false;
  _assemble = true; //must reassemble matrix


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
  
  std::cout << "(ETB) Tight-Binding calculations" << std::endl;

  ModelOptions options;

  if (_assemble && _upt_options.assemble_H) assemble(options);

  if (!_upt_options.assemble_H) create_dummy_H(); 

  if (_upt_solver_options.read_states)
  {
    Messages::info("(ETB) reading old states");

    inst->set_num_states(_upt_solver_options.n_vb, _upt_solver_options.n_cb);

    int n_vb, n_cb;
    inst->read_old_states(_upt_options.load_path, n_vb, n_cb);

    Messages::info("...done\n");
    
    _upt_solver_options.start_vb = n_vb + 1;
    _upt_solver_options.start_cb = n_cb + 1;

    if (_upt_solver_options.n_vb < n_vb) _upt_solver_options.n_vb = n_vb;
    if (_upt_solver_options.n_cb < n_cb) _upt_solver_options.n_cb = n_cb;

    if ( (_upt_solver_options.n_vb > n_vb || _upt_solver_options.n_cb > n_cb) && 
         !_upt_options.assemble_H)
    {
      Messages::warning("Cannot solve more states because assemble_hamiltonian=false");
      _upt_solver_options.n_vb = n_vb;
      _upt_solver_options.n_cb = n_cb;
    }
      
    int num_ev = _upt_solver_options.n_vb + _upt_solver_options.n_cb;
    Complex matel;

    Messages::info("\n(ETB) consistency check:");

    for(int i=1; i<= num_ev; i++)
    {
	matel = inst->get_matel(i,i);
	std::cout << "  <"<<i<<"| H |"<<i<<"> = " << matel << std::endl;
    }
  }

  if (_upt_solver_options.solver.compare("upt_lanczos") == 0) 
  {

    Messages::info("\n(ETB) solving using lanczos");

    if (_upt_solver_options.twice_cb)
      std::cout << "(ETB) twice cb: true " << std::endl;

    if (_upt_solver_options.twice_vb)    
      std::cout << "(ETB) twice vb: true " << std::endl;

    inst->lanczos_diag(_upt_solver_options.start_vb, _upt_solver_options.start_cb, 
		     _upt_solver_options.n_vb, _upt_solver_options.n_cb,
                     _upt_solver_options.guess_vb, _upt_solver_options.guess_cb,
                     _upt_solver_options.min_iter, _upt_solver_options.long_iter,
                     _upt_solver_options.max_iter, _upt_solver_options.fast_tol,
                     _upt_solver_options.long_tol, _upt_solver_options.ort_tol,
                     _upt_solver_options.twice_vb, _upt_solver_options.twice_cb,
		     _upt_solver_options.dynamic);
  }

  if (_upt_solver_options.solver.compare("feast") == 0) 
  {
      Messages::info("Solving Tight Binding with FEAST eigensolver");
      inst->feast(_upt_solver_options.e_min, _upt_solver_options.e_max, _upt_solver_options.m0);
  } 

  if (_upt_solver_options.solver.compare("slepc") == 0) 
  {
      Messages::info("Solving Tight Binding with SLEPc eigensolver");
      //copy_H_to_solver();

  }
    
  Messages::info(" "); 
  Messages::info("(ETB) copy states from uptight"); 

  int hdim = inst->get_H_dim();
  int num_vb = _upt_solver_options.n_vb;
  int num_ev = _upt_solver_options.n_vb + _upt_solver_options.n_cb;


  double *eigvals = new double[num_ev];
  Complex *eigvects = new Complex[hdim*num_ev];
  int *particles = new int[num_ev];
  Complex *eigtmp = eigvects;  // walking pointer on eigenvectors array

  inst->get_states(num_ev,hdim,eigvals,eigvects,particles);

  _solution.resize(num_ev);

  //Set _solution_size in base class TightBinding 
  _solution_size = _solution.size();

  for(int i=0; i< num_ev; i++)
  {

    _solution[i].statistics = "Fermi";
    _solution[i].eigen_energy = eigvals[i];
    _solution[i].temperature = _upt_options.temperature;
    _solution[i].eigen_vector.resize(hdim);

    eigtmp = eigvects + i*hdim;

    for(int j=0; j<hdim;j++)
    {
      _solution[i].eigen_vector[j] = *(eigtmp+j);
    }
     
    if(particles[i]==1)
    { 
       _solution[i].particle = "el";
       if (_upt_options.potential_flag)
       {
	  _solution[i].electro_chem_pot = calculate_fermi_averaged(i);
       }
       else
       {
          _solution[i].electro_chem_pot = _upt_options.el_chem_pot; 	
       }
    }
    else if(particles[i]==-1) 
    {
       _solution[i].particle = "hl";
       if (_upt_options.potential_flag)
       {
	  _solution[i].electro_chem_pot = calculate_fermi_averaged(i);
       }
       else
       {
          _solution[i].electro_chem_pot = _upt_options.hl_chem_pot; 	
       }
    }
    else
    {
       throw SolveFailedException("ETB: unkown particle type");
    }

  }

  delete eigvals;
  delete eigvects;
  delete particles;


  // write state infos on screen.
  write_states();


  //_eigenvector_mag.insert(make_pair(1, std::vector<double>(_N_without_H)));
  for (unsigned int i = 0; i < _solution_size; i++)
  { 
    _eigenvector_mag.insert(make_pair(i, std::vector<double>(_N_without_H)));
    compute_eigenvector_mag(i, _eigenvector_mag[i]); 
  }

  string units("/cm");
  if (_dim == 2) units = "/cm^2";
  else if (_dim == 3) units = "/cm^3";
  declare_solution(MeshStates, NTUPLE, CELL, "1"+units, _solution_size);
  declare_solution(MeshStatesNodes, NTUPLE, NODES, "1"+units, _solution_size);

}

//-------------------------------------------------------------------------
void ETB::do_assemble(const ModelOptions& options)
{


  if( options.get_option("P_matrix",false) )
  {
    int poldir = options.get_option("poldir",0);
    //inst->set_verbose(0);
    //optical matrix computed in UPPER format (now also full works)
    char* sparse_fmt = new char[UPT_MC];    
    memset(sparse_fmt, UPT_PADCHAR, UPT_MC);
    std::string string_fmt = "upper"; 
    string_fmt.copy( sparse_fmt, string_fmt.size() );

    inst->compute_P_matrix(poldir, sparse_fmt);
    //inst->set_verbose(_upt_options.verbose);
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

    inst->compute_H(_upt_options.sparse_fmt);

  }

  _assemble = false;

  _haveS = false;

}

//-------------------------------------------------------------------------
void ETB::print_H(const std::string& outpath) const
{
   char* path; path=new char[UPT_LC];   
   memset(path, UPT_PADCHAR, UPT_LC);
   outpath.copy(path,outpath.size());
   inst->set_workpath(path);
   inst->print_H();
   delete path;
}

//-------------------------------------------------------------------------
void ETB::create_dummy_H(void)
{
   unsigned int nrow;	
   char fmt;

   // get the right value for nrow 
   nrow = compute_H_dim();
   fmt='F';
   // create dummy null matrix with Fortran base 1 indexing
   std::vector<Complex> A(nrow,0.0);
   std::vector<int> JA(nrow,0);
   for(int j=0; j<nrow; j++){ JA[j]=j+1; }
   std::vector<int> IA(nrow+1,0);
   for(int j=0; j<nrow+1; j++){ IA[j]=j+1; }

   inst->set_H_csr(nrow,fmt,A,JA,IA);
 
   // This is needed to create arrays like _el_chem_pot[]
   // for later calculations of averaged Fermi levels
   if (_upt_options.potential_flag)  add_pot_shifts(); 

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
  //if (i_particle == "el" || i_particle == "electron")
  //{
  //  i = i + _upt_solver_options.n_vb;
  //}
  //if (j_particle == "el" || j_particle == "electron")
  //{
  //  j = j + _upt_solver_options.n_vb;
  //}

  //std::cout << "(" << i << "-" << j <<")"<<std::flush;

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
ETB::plot_atomisticdata(void)
{
  //Calculate electron and holes charge density on atoms (hydrogen not included)
  Messages::info("(ETB) compute atom-projected charges"); 

  _el_atomic_charges.resize(_N_without_H, 0.0);
  _hl_atomic_charges.resize(_N_without_H, 0.0);

  compute_atomic_charges("el", _el_atomic_charges);
  compute_atomic_charges("hl", _hl_atomic_charges);

  //Print for debug charges on atoms
  Messages::info("(ETB) print atom-projected charges on files"); 

  string out_path = get_output_directory();
  double* charges;
  charges = new double[get_atomistic_structure()->get_N_atoms()];
  for (unsigned int i = 0; i < _N_without_H; i++) charges[i] = _el_atomic_charges[i];
  for (unsigned int i = _N_without_H + 1; i < get_atomistic_structure()->get_N_atoms(); i++) charges[i] = 0.0;
  get_atomistic_structure()->print_structure(out_path+"/charges_el.xyz", charges);

  for (unsigned int i = 0; i < _N_without_H; i++) charges[i] = _hl_atomic_charges[i];
  for (unsigned int i = _N_without_H + 1; i < get_atomistic_structure()->get_N_atoms(); i++) charges[i] = 0.0;
  get_atomistic_structure()->print_structure(out_path+"/charges_hl.xyz", charges);

  if ( plot_solution("tbstates") || plot_solution("JmolStates") )
  {
    Messages::info("(ETB) write wave functions on files");
    inst->write_states();
  }

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

  std::string sparse_fmt = get_option("sparse_format", "full");
  sparse_fmt.copy(_upt_options.sparse_fmt, sparse_fmt.size() );

  _upt_options.check_bondmap = get_option("check_bondmap", false);

  _upt_options.relat_flag = get_option("relativistic", true);

  _upt_options.temperature = get_option("temperature", SimulationOptions::temperature);

  if (_upt_options.temperature != 0.0)
  {
     Messages::warning("Empirical tight binding calculations are preformed at T=0");
  }

  //_upt_options.opt_flag = get_option("optical_transitions", false);
  //_upt_options.poldir = get_option("polarization_direction", 1);
  _upt_options.opt_flag = false; // these are set via OpticsTB
  _upt_options.poldir = 1;       //   "    "   "   "    "

  //---------------------------------------------------------------------------------------
  // potential_simulation

  _upt_options.potential_sim = get_option("potential_simulation","no_sim");

  if (_upt_options.potential_sim == "no_sim")
  {
   if (!has_option("el_qfermi_level")) 
     Messages::warning("Neither potential_simulation nor el_qfermi_level have been specified");

   if (!has_option("hl_qfermi_level")) 
     Messages::warning("Neither potential_simulation nor hl_qfermi_level have been specified");
  }
  else
  {
    _dd_int = find_simulation(_upt_options.potential_sim);
    if (_dd_int == NULL)
      throw InitFailedException("potential simulation not found");
    _upt_options.potential_flag = true;
  } 

  _upt_options.el_chem_pot = get_option("el_qfermi_level", 0.0);
  _upt_options.hl_chem_pot = get_option("hl_qfermi_level", 0.0);

  //---------------------------------------------------------------------------------------
  // Get strain simulation.
  // both 'strain_model_name' || 'strain_simulation' are accepted keywords
  _upt_options.strain_sim = get_option("strain_model_name", "no_sim");
  _upt_options.strain_sim = get_option("strain_simulation", _upt_options.strain_sim);

  if (_upt_options.strain_sim == "no_sim")
  {
    _upt_options.harrison_flag = get_option("Harrison_scaling", false);
    _upt_options.harrison_flag = get_option("harrison_scaling", _upt_options.harrison_flag);
    _upt_options.d_states_correction = get_option("d_splitting",false); 
  }
  else
  {
    _strain_int.set_simulation(_upt_options.strain_sim);
    _upt_options.harrison_flag = get_option("Harrison_scaling", true);
    _upt_options.d_states_correction = get_option("d_splitting",true);     
  }


  //---------------------------------------------------------------------------------------
  // Dangling bond scaling
  _upt_options.dg_scale = get_option("dangling_bond_scaling",100.0);
  _upt_options.dg_onsite = get_option("dangling_bond_onsite",-200.0);    

  //---------------------------------------------------------------------------------------
  //Choose passivation model
  std::string passivation_model = get_option("passivation_model","hybrid");
  if ( passivation_model == "hybrid" )
      {_upt_options.hybrid_passivation = true;}

  //---------------------------------------------------------------------------------------
 
  _upt_options.assemble_H = get_option("assemble_hamiltonian",true);

  _upt_solver_options.read_states = solopts.get_option("load_states", false);
  std::string load_path = solopts.get_option("load_path", ".");
  load_path.copy(_upt_options.load_path, load_path.size() );
 
  // Solver options: "upt_lanczos"  
  _upt_solver_options.solver = solopts.get_option("solver", "upt_lanczos");
  std::string solver_type = solopts.get_option("solver_type", "cpu");
  if ( solver_type == "cpu") _upt_solver_options.solver_flag = 0;
  if ( solver_type == "gpu") _upt_solver_options.solver_flag = 1;
  if ( solver_type == "gpu-split") _upt_solver_options.solver_flag = 2;

  _upt_solver_options.n_vb =  solopts.get_option("num_valence_eigenvalues", 0);
  if( _upt_solver_options.n_vb == 0) {
    _upt_solver_options.n_vb =  solopts.get_option("num_hole_states", 0);
  }
  
  if (_upt_solver_options.n_vb%2==1) _upt_solver_options.n_vb += 1;
  
  _upt_solver_options.n_cb =  solopts.get_option("num_conduction_eigenvalues", 0);
  if( _upt_solver_options.n_cb == 0) {
    _upt_solver_options.n_cb =  solopts.get_option("num_electron_states", 0);
  }

  if (_upt_solver_options.n_cb%2==1) _upt_solver_options.n_cb += 1;
  
  _upt_solver_options.min_iter =  solopts.get_option("min_iter", 30);
  _upt_solver_options.long_iter =  solopts.get_option("long_iter", 32);
  _upt_solver_options.max_iter =  solopts.get_option("max_iter", 10000);

  bool flag = solopts.get_option("remove_folded_sols_conduction",false);
  if (flag) _upt_solver_options.twice_cb = 1;
  flag = solopts.get_option("remove_folded_sols_valence",false);
  if (flag) _upt_solver_options.twice_vb = 1;
  flag = solopts.get_option("dynamic_search",true);
  if (flag) _upt_solver_options.dynamic = 1;

  //Feast options
  _upt_solver_options.e_min =  solopts.get_option("Emin", 0.0);
  _upt_solver_options.e_max =  solopts.get_option("Emax", 3.0);
  _upt_solver_options.m0 =  solopts.get_option("subspace", 100);

  //---------------------------------------------------------------------------------------
  // output wavevetors format
  std::string out_fmt = get_option("jmol_output_format", "cube");

  // the next lines are ludicrous in order to solve a crazy problem 
  // in passing char* to fortran  
  if (out_fmt=="jvxl") _upt_options.out_format = JVXL;
  else if (out_fmt=="cube") _upt_options.out_format = CUBE;
  else throw InitFailedException("cannot recognize format "+out_fmt);
  
  _upt_options.grid_step = get_option("jmol_grid_step", 0.5);


  //---------------------------------------------------------------------------------------
  
  _vb_shift = 0.0;
  
  //---------------------------------------------------------------------------------------

  _upt_options.band_shift_flag = get_option("add_band_shifts", true);
  _upt_solver_options.guess_vb = solopts.get_option("guess_valence", 0.0);
  _upt_solver_options.guess_cb = solopts.get_option("guess_conduction", 0.0);

  // da togliere e leggere dal database: shift della banda di valenza (che e` 0)
  //_upt_options.vb_shift = options.get_option("vb_shift", 0.0);

  _upt_solver_options.fast_tol =  solopts.get_option("fast_tolerance", 1e-1);
  _upt_solver_options.long_tol =  solopts.get_option("long_tolerance", 1e-9);
  _upt_solver_options.ort_tol =  solopts.get_option("orthogonality_tolerance", 1e-6);

  //Get projection_length for quantum charge projection (Ang)
  _upt_options.projection_length = get_option("projection_length", 2.0);

  //std::cout << "Projection lenght set to " <<  _upt_options.projection_length << std::endl;
  //std::cout << "done" << std::endl;

  // get kpoint
  RealVectorValue k_vec;
  get_option("k_vector",k_vec);
  set_k_point(k_vec);

  // get kpoint as parameter
  Point kp(get_k_point(true));
  get_parameter("k_x", _upt_options.k_point[0], kp(0) );
  get_parameter("k_y", _upt_options.k_point[1], kp(1) );
  get_parameter("k_z", _upt_options.k_point[2], kp(2) );


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
  AtomisticStructure* as = get_atomistic_structure();
  bool randomalloy = as->is_random_alloy();

  unsigned int N = as->get_N_atoms();
  _band_shift.clear();
  _band_shift.resize(N, 0.0);

  const std::vector< Atom >& atom = as->get_structure_atoms();

  const BondMap& b_map = as->get_bond_map();

  if (randomalloy)
  {
    for (unsigned int i = 0; i < N; i++)
    {
      unsigned int nn =  b_map[i].size();
      unsigned int ctr = 0;
      double shift = 0.0;

      for (unsigned int n = 0; n < nn; ++n)
      {
        unsigned int j = b_map[i][n];
        if (atom[j].get_specie() != Specie::H)
        {
          shift += _map_pairs_Evb[
                make_pair(atom[i].get_specie(), atom[j].get_specie())];
          ctr++;
        }
      }
      shift /= ctr;
      _band_shift[i] = -shift + _vb_shift;

    }

    // for the H atoms take the neighbours' values
    for (unsigned int i = 0; i < N; i++)
      if (atom[i].get_specie() == Specie::H)
        _band_shift[i] = _band_shift[b_map[i][0]];
  }
  else
  {
    for (unsigned int i = 0; i < N; i++)
    {

      if (atom[i].get_specie() == Specie::H)
      {
        // the sign is inverted because in upt the potential is subtracted
        _band_shift[i]= - _map_ID_Evb[atom[b_map[i][0]].get_region_ID()] + _vb_shift;
      }
      else
      {
        _band_shift[i]= - _map_ID_Evb[atom[i].get_region_ID()] + _vb_shift;
      }
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
             << _pot_shift[i]  << "  "
             << _el_chem_pot[i] << "  "
             << _hl_chem_pot[i] << std::endl;
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
// Relys on the fact that H are at the end !!
void
ETB::project_atom_strain(void)
{  
  const std::vector<Atom>& structure = get_atomistic_structure()->get_structure_atoms();
  unsigned int Number_of_atoms =  get_atomistic_structure()->get_N_without_H();
  std::vector<double> exx(Number_of_atoms, 0.0);
  std::vector<double> eyy(Number_of_atoms, 0.0);  
  std::vector<double> ezz(Number_of_atoms, 0.0);

  Tensor2Sym epsilon;

  for (unsigned int i = 0; i < Number_of_atoms ; i++)
  { 
    
    if (structure[i].get_elem() != NULL)
    {
      const Elem* el = structure[i].get_elem();

      _strain_int.get_strain(el, el->centroid(), epsilon);

      exx[i] = epsilon(1,1);
      eyy[i] = epsilon(2,2);
      ezz[i] = epsilon(3,3);
    }
  
  }
 
  inst->set_strain(exx,eyy,ezz);
  
}

unsigned int
ETB::compute_H_dim(void)
{
  unsigned int N_atoms_wo_H = get_atomistic_structure()->get_N_without_H();
  unsigned int sum, k_at;

  k_at = 0;
  sum = 0;

  for (unsigned int j = 0; j < N_atoms_wo_H; j++)
  {
    sum += _ion_num_orbitals[j];	  
  }

  return sum;
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
      //cout<<atom_sum<<"  ";
      k_at = k;
      sum +=  _el_chem_pot[j]*atom_sum;
    }
  
  }

  if(_solution[i].particle == "hl" || _solution[i].particle == "hole")
  {

    //cout<<"state hl: "<<i<<endl;    
    for (j = 0; j < N_atoms_wo_H; j++)
    {
      atom_sum = 0.0;
      for (k = k_at; k < k_at + _ion_num_orbitals[j]; k++)
      {
	atom_sum += std::norm(_solution[i].eigen_vector[k]);
      }
      //cout<<atom_sum<<"  ";
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
    
    densatm[j] = atom_sum;

    //cout << j <<"   "<<densatm[j]<<endl;

  }
}

void ETB::get_bulk_edges(void)
{
  
  AtomisticStructure* as = get_atomistic_structure();
  std::set<ID> IDs = as->get_IDset();
  bool randomalloy = as->is_random_alloy();

  for(std::set<ID>::iterator reg = IDs.begin(); reg != IDs.end(); reg++)
  {
      const Material* mat = as->get_device()->get_material( (*reg) );
      {
        if (randomalloy)
        {
          // we get the pure materials' shifts for all atom pairs
          // this works only for binaries up to now
          if (mat->is_alloy())
          {
            const Alloy* alloy = dynamic_cast<const Alloy*>(mat);
            vector<const Material*> mats(2);
            mats[0] = alloy->get_component_A();
            mats[1] = alloy->get_component_B();
            for (int i = 0; i < mats.size(); ++i)
            {
              Database db = mats[i]->get_database();
              db.set_section("atomistic_structure");
              string atom1 = db.get("specie_1", "", true);
              string atom2 = db.get("specie_2", "", true);

              db.set_section("valenceband");
              double vb = db.get("E_v",0.0);

              _map_pairs_Evb[make_pair(Specie(atom1), Specie(atom2))] = vb;
              _map_pairs_Evb[make_pair(Specie(atom2), Specie(atom1))] = vb;
            }
          }
          else
          {
            Database db = mat->get_database();
            db.set_section("atomistic_structure");
            string atom1 = db.get("specie_1", "", true);
            string atom2 = db.get("specie_2", "", true);

            db.set_section("valenceband");
            double vb = db.get("E_v",0.0);

            _map_pairs_Evb[make_pair(Specie(atom1), Specie(atom2))] = vb;
            _map_pairs_Evb[make_pair(Specie(atom2), Specie(atom1))] = vb;
          }
        }


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

double
ETB::get_band_edge(const std::string& edge)
{
  if (edge!="Ec" && edge!="Ev") 
     throw RuntimeException("invalid band-edge"); 

  double sign = (edge=="Ec") ? 1.0 : -1.0;
  double band_edge = sign * numeric_limits<double>::max();

  if (_upt_options.potential_flag)
  {

    Messages::info("(ETB) get band extrema from "+_upt_options.potential_sim);

    ID id = _dd_int->get_solution_id(edge);

    if ((id == INVALID_ID))
      throw RuntimeException("Simulation \'" + _dd_int->get_name() +
          "\' lacks band edge solution variables.");

    vector<double> edges(8,0.0);

    MeshBase::const_element_iterator it = get_mesh().active_elements_begin();
    const MeshBase::const_element_iterator end = get_mesh().active_elements_end();

    for ( ; it != end; ++it)
    { 
      const Elem* elem = *it;
      vector<Point> p(elem->n_nodes());

      for (size_t i = 0; i < elem->n_nodes(); ++i)
        p[i] = elem->point(i);

      _dd_int->get_solution(elem, id, edges, p);

      for (size_t i = 0; i < elem->n_nodes(); ++i)
      {
        double be = edges[i];
        // sign trick: 
        // for Ec::  be<edge ? be:edge
        // for Ev:: -be<-edge => be>edge ? be:edge
        band_edge = (sign*be < sign*band_edge) ? be : band_edge;
      }

    }

  }
  else
  {
    //computes guess from database band edges
    std::set<ID> IDs = get_atomistic_structure()->get_IDset();
    std::set<ID>::iterator reg;

    for (reg = IDs.begin(); reg != IDs.end(); reg++)
    {
      if (edge=="Ec")
      {
        band_edge = (_map_ID_Ecb[*reg] < band_edge) ? _map_ID_Ecb[*reg] : band_edge; 
      }
      else
      {
        band_edge = (_map_ID_Evb[*reg] > band_edge) ? _map_ID_Evb[*reg] : band_edge; 
      }
    }
  }

  return band_edge;

}


void
ETB::get_band_extrema(double& cb_min, double& vb_max)
{
  if (_upt_options.potential_flag)
  {
    Messages::info("(ETB) get band extrema from "+_upt_options.potential_sim);

    string cbedge("Ec");
    string vbedge("Ev");

    // here we assume that the simulation exists (it was checked before)
    ID cb_id = _dd_int->get_solution_id(cbedge);
    ID vb_id = _dd_int->get_solution_id(vbedge);

    if ((cb_id == INVALID_ID) || (vb_id == INVALID_ID))
      throw RuntimeException("Simulation \'" + _dd_int->get_name() +
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

      _dd_int->get_solution(elem, vb_id, vb_edges, p);
      _dd_int->get_solution(elem, cb_id, cb_edges, p);

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

  //CELL values
  //
  if (values.count(ElQuantumDensity))
    {
    if (_dim == 3)
      values[ElQuantumDensity][0] = build_rho3d(_el_atomic_charges, elem, elem->centroid());
    else if (_dim == 2)
      values[ElQuantumDensity][0] = build_rho2d(_el_atomic_charges, elem, elem->centroid());
    else if (_dim == 1)
      values[ElQuantumDensity][0] = build_average_rho1d(_el_atomic_charges, elem);
    }

  if (values.count(HlQuantumDensity))
    {
    if (_dim == 3)
      values[HlQuantumDensity][0] = build_rho3d(_hl_atomic_charges, elem, elem->centroid());
    else if (_dim == 2)
      values[HlQuantumDensity][0] = build_rho2d(_hl_atomic_charges, elem, elem->centroid());
    else if (_dim == 1)
      values[HlQuantumDensity][0] = build_average_rho1d(_hl_atomic_charges, elem);
    }

  if (values.count(MeshStates))
    {
    for (unsigned int i = 0; i < _solution_size; i++)
      {
      if (_dim == 3)
        values[MeshStates][i] = build_rho3d(_eigenvector_mag[i], elem, elem->centroid());
      else if (_dim == 2)
        values[MeshStates][i] = build_rho2d(_eigenvector_mag[i], elem, elem->centroid());
      else if (_dim == 1)
        values[MeshStates][i] = build_average_rho1d(_eigenvector_mag[i], elem);
      }
    }

  //NODES values
  unsigned int np = p.size();


  //Get point coordinate as physical_point
  Point phys_p;

  if (values.count(ElQuantumDensityNodes))
    {
    for (unsigned int n = 0; n < np; n++)
      {
      if (_dim == 3)
        {
        phys_p = FE< 3, libMeshEnums::LAGRANGE>::map(elem, p[n]);
        values[ElQuantumDensityNodes][n] = build_rho3d(_el_atomic_charges, elem, phys_p);
        }
      else if (_dim == 2)
        {
        phys_p = FE< 2, libMeshEnums::LAGRANGE>::map(elem, p[n]);
        values[ElQuantumDensityNodes][n] = build_rho2d(_el_atomic_charges, elem, phys_p);
        }
      else if (_dim == 1)
        {
        throw ModelErrorException("ElQuantumNodes is not supported in 1D "
            "calculations. Use ElQuantumDensity (CELL variable)");
        }
      }
    }

  if (values.count(HlQuantumDensityNodes))
    {
    for (unsigned int n = 0; n < np; n++)
      {
      if (_dim == 3)
        {
        phys_p = FE< 3, libMeshEnums::LAGRANGE>::map(elem, p[n]);
        values[HlQuantumDensityNodes][n] = build_rho3d(_hl_atomic_charges, elem, phys_p);
        }
      else if (_dim == 2)
        {
        phys_p = FE< 2, libMeshEnums::LAGRANGE>::map(elem, p[n]);
        values[HlQuantumDensityNodes][n] = build_rho2d(_hl_atomic_charges, elem, phys_p);
        }
      else if (_dim == 1)
        {
        throw ModelErrorException("HlQuantumNodes is not supported in 1D "
            "calculations. Use HlQuantumDensity (CELL variable)");
        }

      }
    }

  if (values.count(MeshStatesNodes))
    {
    for (unsigned int n = 0; n < np; n++)
      {

      for (unsigned int i = 0; i < _solution_size; i++)
        {
        if (_dim == 3)
          {
          phys_p = FE< 3, libMeshEnums::LAGRANGE>::map(elem, p[n]);
          values[MeshStatesNodes][_solution_size * n + i] = 
            build_rho3d(_eigenvector_mag[i], elem, phys_p);
          }
        else if (_dim == 2)
          {
          phys_p = FE< 2, libMeshEnums::LAGRANGE>::map(elem, p[n]);
          values[MeshStates][_solution_size * n + i] = 
            build_rho2d(_eigenvector_mag[i], elem, phys_p);
          }
        else if (_dim == 1)
          {
          throw ModelErrorException("MeshStatesNodes is not supported in 1D "
              "calculations. Use MeshStates (CELL variable)");
          }

        }
      }
    }

}


double
ETB::build_rho3d(const std::vector<double>& tb_density, const Elem* elem, const Point& r)
{
  double scale = get_atomistic_structure()->get_scale();
  //double tau = 1.0 / _upt_options.projection_length; // projection in Angstroms
  const double sigma = _upt_options.projection_length;
  const double sigma2 = 2.0*sigma*sigma;

  double deltar, deltar2;
  double rho = 0.0;
  double x1, y1, z1;
  double x ,y, z;
  
  x = r(0)*scale; y = r(1)*scale; z = r(2)*scale;

  if (tb_density.size() != _N_without_H)
  {
    Messages::error("ERROR IN ETB: trying to build mesh density but tb density has wrong size");
    exit(1);
  }
 
  //! get structure
  const std::vector<Atom>& structure = get_atomistic_structure()->get_structure_atoms();
  //! get atoms that contribute to the density of an element
  //std::cout<<"(ETB) get neigh_atoms "<<endl;
  const std::vector<unsigned int>& atoms = get_elem_atoms(elem->id());

  //std::cout<<"point: "<<x<<" "<<y<<" "<<z<<endl;

  for (unsigned int id = 0; id  < atoms.size(); id++)
  {
    unsigned int iatm = atoms[id];

    //cout<<iatm<<":";
    
    // position is in Angstrom
    x1 = structure[iatm].get_position(0);
    y1 = structure[iatm].get_position(1);
    z1 = structure[iatm].get_position(2);
    
    //delta_r is in Angstrom
    //deltar = sqrt( (x - x1) * (x - x1) + (y - y1) * (y - y1) + (z - z1) * (z - z1));
    deltar2 = (x - x1) * (x - x1) + (y - y1) * (y - y1) + (z - z1) * (z - z1);
    //cout<<deltar<<";";

    //rho = rho + (tb_density[iatm] * exp(-deltar * tau));      
    rho = rho + (tb_density[iatm] * exp(- deltar2/sigma2 )); 

  }
  
  //rho = rho *tau*tau*tau / (8.0 * 3.141592653589793);

  double norm = 2.0 * 3.141592653589793 * sigma;
  
  rho = rho / sqrt(norm*norm*norm);

  //scale rho from q/Ang^3 to q/cm^3 (carrier density)
  rho =  rho * 1e24;
  return rho;

}


double
ETB::build_rho2d(const std::vector<double>& tb_density, const Elem* elem, const Point& r)
{ 
  double scale = get_atomistic_structure()->get_scale();
  //double tau = 1.0 / _upt_options.projection_length;
  double deltar, deltar2, uhatom;
  const double sigma = _upt_options.projection_length;
  const double sigma2 = 2.0*sigma*sigma;  
  double rho = 0.0;
  double x1, y1, z1;
  double x ,y, z;


  x = r(0)*scale; y = r(1)*scale; z = 0.0;

  if (tb_density.size() != _N_without_H)
  {
    std::cerr << "ERROR IN ETB: trying to mesh density "
      "but tb_density has wrong size" << std::endl;
    exit(1);
  }

  const std::vector<Atom>& structure = get_atomistic_structure()->get_structure_atoms();

  //! get atoms that contribute to the density of an element
  const std::vector<unsigned int>& atoms = get_elem_atoms(elem->id());


  for (unsigned int id = 0; id  < atoms.size(); id++)
  {
    unsigned int iatm = atoms[id];    

    x1 = structure[iatm].get_position(0);
    y1 = structure[iatm].get_position(1);
    z1 = 0.0;
    
    deltar2 = (x - x1) * (x - x1) + (y - y1) * (y - y1);
    
    //rho = rho + (tb_density[iatm] * tau * tau * exp(-1.0 * deltar * tau));
    rho = rho + (tb_density[iatm] * exp(- deltar2/sigma2 ));

  }

  double norm = 2.0 * 3.141592653589793 * sigma; 
  rho = rho / norm;

  //scale rho from q/Ang^2 to q/cm^2 (carrier density)
  rho =  rho * 1e16;

  return rho;

}


double
ETB::build_average_rho1d(const std::vector<double>& tb_density, const Elem* elem)
{
  //Projection length in mesh units
  double scale = get_atomistic_structure()->get_scale();
  double tau = 1.0 / _upt_options.projection_length;
  const double deltar_max = 10 / tau; //Maximum cutoff distance
  double deltar, uhatom, l;
  double rho = 0.0;
  double x_atm, x1, x2;

  assert(elem->n_nodes() == 2);

  x1 = elem->point(0)(0)*scale; x2 = elem->point(1)(0)*scale;

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

  const std::vector<Atom>& structure = get_atomistic_structure()->get_structure_atoms();

  for (unsigned int iatm = 0; iatm  < _N_without_H; iatm++)
    {

      x_atm = structure[iatm].get_position(0);

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

  //scale rho from q/Ang to q/cm (carrier density)
  rho =  rho * 1e8;

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
  declare_solution(ElQuantumDensityNodes, REAL, NODES, "q"+units);
  declare_solution(HlQuantumDensityNodes, REAL, NODES, "q"+units);
  declare_solution(MeshStates, NTUPLE, CELL, "1"+units, 1);
  declare_solution(MeshStatesNodes, NTUPLE, NODES, "1"+units, 1);

  if (plot_solution("ProbabilityDensity"))
	    add_plot_variable(MeshStates);

  if (plot_solution("ProbabilityDensityNodes"))
            add_plot_variable(MeshStatesNodes);

}

// -------------------------------------------------------------------------------------------------

void ETB::do_copy_H_to_solver( )
{

  int size_matrix = inst->get_H_dim();
  
  int non_zeros_number[size_matrix];

  for (int row = 0 ; row < size_matrix; row++)	
  {
      non_zeros_number[row] = inst->get_H_row_size(row);	  
  }

  EigenSolver::preallocate_H_matrix(size_matrix,  non_zeros_number);


  //----------------------------------------------------------------------------------------------------//
 
  //write data of columns in each row
  
  for (int row = 0 ; row < size_matrix; row++)
  {
      Complex *row_vals;
      int *cols;
      int n_cols = 0;
	
      vector<unsigned int> column_vector;
      vector<Complex> row_values;

      n_cols = inst->get_H_row_size(row);
 
      cols = new int[n_cols]; //PetscInt[n_cols];
      row_vals = new Complex[n_cols];

      inst->get_H_row(row, cols, row_vals);

      column_vector.resize(n_cols);
      column_vector.clear();

      row_values.resize(n_cols);
      row_values.clear();

      for (int i = 0; i < n_cols; i++)
      {
	  column_vector[i] = cols[i];
	  row_values[i] = row_vals[i];
      }
     
      EigenSolver::insert_H_row( row, column_vector, row_values);
  }

  EigenSolver::finalize_H_assembly();
  
}

int
ETB::get_H_dim() const
{
  return inst->get_H_dim();
}

int
ETB::get_H_nnz() const
{
  return inst->get_H_nnz();
}

void
ETB::get_H_csr(std::vector<Complex>& A,std::vector<int>& JA, std::vector<int>& IA) const
{
   int nrow = inst->get_H_dim();
   inst->get_H_csr(nrow,' ',A,JA,IA);
}

unsigned int
ETB::get_number_of_bands(void) const
{
   unsigned int num = _upt_options.relat_flag ? 20 : 10;
   return num;
}

void
ETB::setup_atomistic_structure(void)
{
  string name(get_option("atomistic_structure", ""));
  if (!name.empty())
  {
    _atomistic_structure = get_environment().get_device().get_atomistic_structure(name);
    if (_atomistic_structure == NULL)
      throw InitFailedException("No atomistic structure \'" + name + "\' found "
          "for simulation \'" + get_name());

    if (!includes_regions(_atomistic_structure->get_IDset())) 
    { 
       Messages::warning("Module will restrict the atomistic structure");      


       AtomisticStructure* new_as = AtomisticStructure::create(*get_atomistic_structure());
       

       _atomistic_structure = new_as;
    
    }
  }
}


