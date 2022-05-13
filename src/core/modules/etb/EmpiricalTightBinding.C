// $Id$

#include <boost/filesystem/operations.hpp>

#include "EmpiricalTightBinding.h"
#include "PhysicalModel.h"
#include "AtomisticStructure.h"
#include "Atom.h"
#include "SimulationOptions.h"
#include "TiberLinearSystem.h"
#include "TiberCad.h"
#include "UptWrapper.h"
#include "uptight.h"
#include "Material.h"
#include "Alloy.h"
#include "Messages.h"
#include "EigenSolver.h"
#include "RotatedCrystal.h"
#include "Utils.h"

#include "libmesh/mesh.h"
#include "libmesh/dof_map.h"

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
using namespace libMesh;


namespace
{
  bool compare_eigen_energy(const EigenvalueProblem::eigen_state& state1,
      const EigenvalueProblem::eigen_state& state2)
  {
    if ((state1.particle == "hl") && (state1.particle == state2.particle))
      return (state2.energy < state1.energy);

    return(state1.energy < state2.energy);
  }
}


//--------------------------------------------------------------

ETB::ETB(const ModelOptions& options)
: TightBinding(options),
  _upt_options()
{
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
 pdos_flag(false),
 potential_flag(0),
 opt_flag(0),
 poldir(1),
 hybrid_passivation(false),
 dg_scale(1.0),
 dg_onsite(-200.0),
 grid_step(0.5),
 compute_densities(false)
{
  //c_axis.reserve(3);
  c_axis[0]=0.0; c_axis[1]=0.0; c_axis[2]=1.0;
  //k_point.reserve(3);
  k_point[0]=0.0; k_point[1]=0.0; k_point[2]=0.0;
  default_path = new char[UPT_LC];  memset(default_path, UPT_PADCHAR, UPT_LC);
  database_path = new char[UPT_LC]; memset(database_path, UPT_PADCHAR, UPT_LC);
  work_path = new char[UPT_LC];     memset(work_path, UPT_PADCHAR, UPT_LC);
  load_path = new char[UPT_LC];     memset(load_path, UPT_PADCHAR, UPT_LC);
  write_state = new char[UPT_LC];   memset(write_state, UPT_PADCHAR, UPT_LC);
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
  delete[] write_state;
  delete[] out_path;
  delete[] upt_filename;
  delete[] gen_outfile;
  delete[] sparse_fmt;
  //delete[] out_format;
  delete[] database_path;
  delete[] default_path;
}

ETB::UptSolverOptions::UptSolverOptions(void)
  :solver("lanczos"),
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
   dynamic(0)	
{
}

ETB::UptSolverOptions::~UptSolverOptions(void)
{

}






//-------------------------------------------------------------------------
void
ETB::do_init(void){

  std::cerr << "("+get_name()+") Empirical TB Initialisation..." << std::endl;

  inst = UptWrapper::create();

  inst->set_mpi_comm(this->get_solver_communicator().get());

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
  std::string work_path = get_output_directory();
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

  Messages::info("("+get_name()+") creating map elem->atoms");

  build_map_elem_atoms(_upt_options.projection_length);

  cerr << "done\n";

  Messages::info("("+get_name()+") database path: "+database_path);
  Messages::info("("+get_name()+") default  path: "+default_path);
  Messages::info("("+get_name()+") work path: "+work_path);
  Messages::info("("+get_name()+") output path: "+out_path);
  

  _init = true;      // initialization must be called
  _assemble = true;  // matrix assemble must be done



  // We add a second system just to contain the density
  create_equation_system("linear");
  TiberLinearSystem& linsys = get_equation_system<TiberLinearSystem>(0);
  linsys.add_variable("edens", libMeshEnums::FIRST, &(this->get_region_ids()));
  linsys.add_variable("hdens", libMeshEnums::FIRST, &(this->get_region_ids()));
  linsys.init();

}

//-------------------------------------------------------------------------
void ETB::do_reinit(void)
{

  inst->cleanuptight();

  if (!includes_regions(get_atomistic_structure()->get_IDset()))
  {
     std::cerr << "("+get_name()+") restrict on active regions" << std::endl;

     std::set<ID> reg_ids;
     get_region_ids(reg_ids);

     get_atomistic_structure()->dorestrict(reg_ids);
    
     std::cerr<< "("+get_name()+") Build Map Elem->Atom "<<std::endl;

     double Rmax = build_map_elem_atoms(_upt_options.projection_length);
     
     std::cerr<< "("+get_name()+") Rmax= "<<Rmax<<std::endl;
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

  Point kp(get_k_point(true));

  // create default name for state files
  string statefile("states.upt");
  if (kp.norm() > 1e-6) // k = Gamma
  {
    ostringstream os;
    os.precision(4);
    os << "states_(" << fixed << kp(0) << "," << kp(1) << "," << kp(2) << ").upt";
    statefile = os.str();
  }
  
  ModelOptions& solopts = get_solver_options();
  // set writing state file  
  string write_state = solopts.get_option("states_file", statefile);
  if (write_state.at(0) != '/')
    write_state =  get_output_directory() +  "/" + write_state;
  memset(_upt_options.write_state, UPT_PADCHAR, UPT_LC);
  write_state.copy(_upt_options.write_state, write_state.size());
  
  // set load state file  
  _upt_solver_options.read_states = false;
  if (solopts.find_option("load_file") || solopts.find_option("load_path") ||
      solopts.get_option("load_states", false))
  {
    _upt_solver_options.read_states = true;
    string load_path = solopts.get_option("load_path", get_output_directory());
    string load_file(load_path + "/" + statefile);
  
    load_file = solopts.get_option("load_file", load_file);

    memset(_upt_options.load_path, UPT_PADCHAR, UPT_LC);
    load_file.copy(_upt_options.load_path, load_file.size());
  }

  // Get the minimum CB and maximum VB edges
  if (!solopts.find_option("guess_conduction") ||
      !solopts.find_option("guess_valence"))
  {
    double cb_min, vb_max;
    get_band_extrema(cb_min, vb_max);
    ostringstream os;
    os << "      CB min = " << cb_min << "  VB max = " << vb_max << endl;
    Messages::info(os.str());

    double gap = cb_min - vb_max;
    if (gap <= 0.0)
    {
      Messages::warning("Your system apparently does not have a global gap: "
          "cannot find reasonable guess.");
      Messages::warning("Will use mean value of band edges");
      cb_min = (cb_min + vb_max) / 2.0;
      vb_max = cb_min;
    }
    else
    {
      cb_min -= gap / 4.0;
      vb_max += gap / 4.0;
    }

    _upt_solver_options.guess_cb = cb_min;
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

    upt_filename = get_output_directory() + "/" +
        get_atomistic_structure()->get_name() + ".upg";

    //Messages::info("("+get_name()+") printing structure "+upt_filename);

    if (get_communicator().rank() == 0)
      get_atomistic_structure()->print_upg(upt_filename, _upt_options.etb_dataset,
                                                      !_upt_options.band_shift_flag);
    get_communicator().barrier();


    std::cout << "("+get_name()+") Number of atoms: " <<get_atomistic_structure()->get_N_atoms() << std::endl;

    std::cout << "("+get_name()+") Numb. without H: " <<get_atomistic_structure()->get_N_without_H()
	      << std::endl;
  }

  std::size_t length = 0;
  length = upt_filename.copy(_upt_options.upt_filename, upt_filename.size() );


  // Use rotated crystal to obtain direction of c-axis (z-direction)
  get_c_axis();

  std::cout << "("+get_name()+") c-axis: "<<_upt_options.c_axis[0]<<" "
            <<_upt_options.c_axis[1]<<" "<<_upt_options.c_axis[2]<<std::endl;  
  

  //  Set parameters for Uptight instance
  inst->set_paths(_upt_options.default_path, _upt_options.database_path, 
                  _upt_options.work_path, _upt_options.out_path);
  
  std::cout << "("+get_name()+") fill parameter " << std::endl;

  inst->fill_param(_upt_options.verbose, _upt_options.upt_filename,
		   _upt_options.gen_outfile, _upt_options.sparse_fmt,
		   _upt_options.max_TB_order, _upt_options.harrison_flag,
		   _upt_options.relat_flag, _upt_options.potential_flag,
		   _upt_options.opt_flag, _upt_options.poldir,
		   _upt_options.c_axis, _upt_options.check_bondmap,
       _upt_options.dg_scale, _upt_options.dg_onsite,
       _upt_options.hybrid_passivation);
 
  std::cout << "("+get_name()+") set solver flag "<< _upt_solver_options.solver_flag << std::endl;
  inst->set_solver_flag(_upt_solver_options.solver_flag); 
	  
  inst->set_output((int) _upt_options.out_format, _upt_options.grid_step);

  //std::cout << "("+get_name()+") fill parameter done" << std::endl;
  //std::cout.flush();

  if(inst->inituptight() != 0){ 
    throw InitFailedException("internal handlers do not match"); }

  _ion_num_orbitals.resize(get_atomistic_structure()->get_N_atoms(), 0);

  inst->get_ion_numorbitals(_ion_num_orbitals);

  _N_without_H = get_atomistic_structure()->get_N_without_H();

  if (has_new_k())
  {
    Point kp(get_k_point(true));
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
    // NOTE: the z axis is assumed to be the reference axis as used in uptight
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




void ETB::do_calculate_density_at_k(DofField& density)
{
  _el_atomic_charges.clear();
  _el_atomic_charges.resize(_N_without_H, 0.0);
  _hl_atomic_charges.clear();
  _hl_atomic_charges.resize(_N_without_H, 0.0);

  compute_atomic_charges("el", _el_atomic_charges);
  compute_atomic_charges("hl", _hl_atomic_charges);

  //cerr << "####\n";
  //for (int i = 0; i < _el_atomic_charges.size(); ++i)
  //{
  //  cerr << i+1 << " " << _el_atomic_charges[i] << " " << _hl_atomic_charges[i] << endl;
  //}
  //cerr << "####\n";

  density = _el_atomic_charges;
  density.reserve(2 * _el_atomic_charges.size());
  density.insert(density.end(), _hl_atomic_charges.begin(), _hl_atomic_charges.end());
}




void ETB::do_solve_for_kpoint(const Point& k_point)
{
  reinit();
  call_uptight();
}


//-------------------------------------------------------------------------
void ETB::do_solve(void) {

  if (_upt_options.compute_densities)
  {
    // calculate the atomic density
    DofField dens;
    integrate_density(dens);

    // put densities back into the qdensity vectors
    unsigned int n_atoms = dens.size() / 2;
    for (unsigned int i = 0; i < n_atoms; ++i)
    {
      //cerr << dens[i] << " " << dens[i + n_atoms] << std::endl;
      _el_atomic_charges[i] = dens[i];
      _hl_atomic_charges[i] = dens[i + n_atoms];
    }

    // The qdens_sys system contains the nodal quantum density
    TiberLinearSystem& qdens_sys = get_equation_system<TiberLinearSystem>();
    NumericVector<Number>& qdens = *qdens_sys.solution;

    FEType fe_type = qdens_sys.variable_type(0);
    UniquePtr<FEBase> fe(build_finite_element(_dim, fe_type));

    DofMap& dof_map = qdens_sys.get_dof_map();
    vector<unsigned int> dof_indices_e;
    vector<unsigned int> dof_indices_h;

    // we play a trick: only DoFs with value == -1 need to be updated
    qdens.zero();
    qdens.add(-1.0);

    // the cutoff distance for near atoms in A
    double cutoff = _upt_options.projection_length;

    Utils::Progress prog("Project densities", get_mesh().n_active_elem());


    MeshBase::const_element_iterator el = this->active_local_elements_begin();
    const MeshBase::const_element_iterator end_el = this->active_local_elements_end();
    for ( ; el != end_el; ++el)
    {
      const Elem* elem = *el;
      dof_map.dof_indices(elem, dof_indices_e, 0);
      dof_map.dof_indices(elem, dof_indices_h, 1);

      for (unsigned int n = 0; n < elem->n_nodes(); n++)
      {
        //if (qdens(dof_indices_e[n]) < 0)
        {
          auto densities = project_densities(elem, elem->point(n), cutoff);
          qdens.set(dof_indices_e[n], densities.first);
          qdens.set(dof_indices_h[n], densities.second);
        }
      }
      prog.progress_message();

    }

    qdens.close();
  }
  else
    call_uptight();
}


void
ETB::call_uptight(void)
{

  // check unused tags in solver_options
  const ModelOptions& sol_opt = get_solver_options();
  //sol_opt.find_option("simulation"); // remove simulation name
  //sol_opt.check_unused();
  
  //std::cout << "(ETB) Tight-Binding calculations" << std::endl;
  
  inst->set_statefile(_upt_options.write_state);

  ModelOptions options;

  if (_assemble && _upt_options.assemble_H) assemble(options);

  if (!_upt_options.assemble_H) create_dummy_H(); 

  if (_upt_solver_options.read_states)
  {
    Messages::info("("+get_name()+") reading old states");

    inst->set_num_states(_upt_solver_options.n_vb, _upt_solver_options.n_cb);

    int n_vb = 0, n_cb = 0;

    // the following is stupid, but needed for the way the path is stored
    string tmpstr;
    for (int i = 0; i < UPT_LC; ++i)
    {
      if (_upt_options.load_path[i] != UPT_PADCHAR)
        tmpstr.append(1, _upt_options.load_path[i]);
      else
        break;
    }

    boost::filesystem::path loadfile(tmpstr);
    if (!boost::filesystem::exists(loadfile))
    {
      Messages::warning("File " + loadfile.string() +
          " does not exist: skipping file reading.");
    }
    else
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

    Messages::info("\n("+get_name()+") consistency check:");

    for(int i=1; i<= num_ev; i++)
    {
	Complex matel = inst->get_matel(i,i);
	std::cout << "  <"<<i<<"| H |"<<i<<"> = " << matel << std::endl;
    }
  }

  string solver_type = get_solver_options().get_option("solver_type", "cpu");
  if (solver_type.compare("slepc") == 0)
  {
    Messages::info("Solving Tight Binding with SLEPc eigensolver");
    initialize_solution_container(_upt_solver_options.n_vb + _upt_solver_options.n_cb);
    solve_eigenvalue_problem(_upt_solver_options.n_vb + _upt_solver_options.n_cb, _upt_solver_options.guess_cb);

    int num_states = _upt_solver_options.n_vb + _upt_solver_options.n_cb;

    for (int i = 0; i < num_states; ++i)
    {
      // +1 is electron, -1 is hole
      int particle = 1;
      if (_solution[i].particle == "hl")
        particle = -1;

      inst->set_state(num_states, i, _solution[i].eigen_vector.size(),
                                     _solution[i].eigen_energy,
                                     _solution[i].eigen_vector, particle);
    }
  }
  else
  {
    if ((_upt_solver_options.solver.compare("upt_lanczos") == 0) ||
        (_upt_solver_options.solver.compare("lanczos") == 0))
    {

      Messages::info("\n("+get_name()+") solving using lanczos");

      inst->lanczos_diag(_upt_solver_options.start_vb, _upt_solver_options.start_cb,
          _upt_solver_options.n_vb, _upt_solver_options.n_cb,
          _upt_solver_options.guess_vb, _upt_solver_options.guess_cb,
          _upt_solver_options.min_iter, _upt_solver_options.long_iter,
          _upt_solver_options.max_iter, _upt_solver_options.fast_tol,
          _upt_solver_options.long_tol, _upt_solver_options.ort_tol,
          _upt_solver_options.dynamic, _upt_solver_options.bitoff);
    }
    else if (_upt_solver_options.solver.compare("feast") == 0)
    {
      Messages::info("Solving Tight Binding with FEAST eigensolver");
      inst->feast(_upt_solver_options.e_min, _upt_solver_options.e_max, _upt_solver_options.m0);
    }
    else if (_upt_solver_options.solver.compare("jd") == 0)
    {
      Messages::info("\n("+get_name()+") solving using Jacobi-Davidson");

      inst->jacobidavidson(_upt_solver_options.start_vb, _upt_solver_options.start_cb,
          _upt_solver_options.n_vb, _upt_solver_options.n_cb,
          _upt_solver_options.guess_vb, _upt_solver_options.guess_cb,
          _upt_solver_options.long_tol);

    }
    else if (_upt_solver_options.solver.compare("lapack") == 0)
    {
      Messages::info("\n("+get_name()+") solving using LAPACK");
      inst->lapack(_upt_solver_options.n_vb, _upt_solver_options.n_cb,
          _upt_solver_options.guess_vb, _upt_solver_options.guess_cb);
    }


    Messages::info(" ");
    Messages::info("("+get_name()+") copy states from uptight");

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

    delete [] eigvals;
    delete [] eigvects;
    delete [] particles;
  }


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


std::pair<unsigned int, double>
ETB::read_slepc_solution(void)
{
  //--------------------------------------------------------------------
  //how many solutions do we have from SLEPC?
  unsigned int number_of_converged_solutions;

  number_of_converged_solutions = EigenSolver::number_of_converged_eigenvalues();

  unsigned int number_of_eigenstates = _upt_solver_options.n_vb + _upt_solver_options.n_cb;
  //--------------------------------------------------------------------
  //read eigenvalues
  //store also eigenvalue index for sorting

  vector<EigenvalueProblem::eigen_state>  ev(number_of_converged_solutions);
  double shift = EigenSolver::get_shift();

  // if we have already solutions, we should use their energy levels to discriminate
  // between el and hl
  if (_upt_solver_options.n_vb < _solution.size() &&
      _solution[_upt_solver_options.n_vb].eigen_vector.size() > 0)
  {
    shift = _solution[_upt_solver_options.n_vb].eigen_energy;
  }
  else if ((_upt_solver_options.n_vb > 0) &&
           (_solution[_upt_solver_options.n_vb - 1].eigen_vector.size() > 0))
  {
    shift = _solution[_upt_solver_options.n_vb - 1].eigen_energy;
  }

  for (unsigned ind = 0; ind < number_of_converged_solutions; ind++)
  {
    ev[ind].energy =  EigenSolver::get_eigenvalue(ind);
    ev[ind].index = ind;

    if (ev[ind].energy > shift)
      ev[ind].particle = "el";
    else
      ev[ind].particle = "hl";
  }

  // sorting of the solutions
  // we sort both electrons and holes by distance from the ground state

  sort(ev.begin(), ev.end(), compare_eigen_energy);

  if (verbose() > 1)
  {
    Messages m;
    ostringstream os;
    os << "converged eigenenergies (" << number_of_converged_solutions
        << "):";
    m.info(os.str());
    m.indent();

    os.str("");
    for (unsigned int i = 0; i < number_of_converged_solutions; ++i)
    {
      os << ev[i].energy << " ";
      if (i%8 == 7)
        os << "\n";
    }
    Messages::info(os.str());
    m.newline();
  }


  // find the first electron state

  unsigned int first_el_index = number_of_converged_solutions;
  bool finish = false;

  for (unsigned int i = 0; i < number_of_converged_solutions; i++)
  {
    if (ev[i].particle == "el")
    {
      first_el_index = i;
      break;
    }
  }
  //--------------------------------------------------------------------
  //read eigenvectors

  // The idea is that arriving here the _solution structure is set up,
  // but all eigenvectors are empty or already calculated, valid eigenstates


  // the first num_hl_states in _solution are holes, the upper
  // num_el_states ones are electrons


  const int n_states = _solution.size();

  //----------------------------------------------------------------------
  for (unsigned int i = 0; i < number_of_converged_solutions; i++)
  {
    // calculate the index in the final solution structure
    int index = static_cast<int>(_upt_solver_options.n_vb + i) - first_el_index;
    if (ev[i].particle == "hl")
      index = i;

    if (((ev[i].particle == "hl") && (index >= _upt_solver_options.n_vb)) ||
        ((ev[i].particle == "el") && (index < _upt_solver_options.n_vb)) ||
        ((index < 0) || (index >= n_states)))
        continue;


    // we need a small delta to decide if two states may be degenerate
    // TODO adjust it automatically
    const double delta = 1e-5;

    // look for the first available slot
    if (ev[i].particle == "el")
    {
      while ((index < n_states) &&
             (_solution[index].eigen_vector.size() > 0))
      {
        index++;
      }

      if (index > _upt_solver_options.n_vb)
      {
        if ((index >= n_states) ||
            (ev[i].energy < (_solution[index - 1].eigen_energy - delta))) // go to the next state
          continue;

      }
    }
    else // if (ev[i].particle == "hl")
    {
      while ((index < static_cast<int>(_upt_solver_options.n_vb)) &&
             (_solution[index].eigen_vector.size() > 0))
      {
        index++;
      }

      if (index > 0)
      {
        if ((index >= static_cast<int>(_upt_solver_options.n_vb)) ||
            (ev[i].energy > _solution[index - 1].eigen_energy + delta)) // go to the next state
          continue;

      }
    }

    // we found a (potentially) valid slot and fill it
    _solution[index].eigen_energy = ev[i].energy;
    _solution[index].particle = ev[i].particle;
    _solution[index].statistics = "Fermi";
    _solution[index].temperature = _upt_options.temperature;


    unsigned int solution_number = ev[i].index;

    //EigenSolver::get_eigen_vector(solution_number, temp);
    EigenSolver::get_eigen_vector(solution_number, _solution[index].eigen_vector);

    this->get_solver_communicator().allgather(_solution[index].eigen_vector);

    if (_upt_options.potential_flag)
    {
      _solution[index].electro_chem_pot = calculate_fermi_averaged(index);
    }
    else
    {
      if (_solution[index].particle == "el")
        _solution[i].electro_chem_pot = _upt_options.el_chem_pot;
      else
        _solution[i].electro_chem_pot = _upt_options.hl_chem_pot;
    }

    //
    //normalization
    //
    //double norm = eigenstate_norm(index);

    //for (unsigned int j = 0; j < number_of_all_dofs; j++)
    //  _solution[index].eigen_vector[j] /= Complex(norm, 0.0);



  }

  // the 1e-5 below is to not make the Hamiltonian singular,
  // and to be sure to take all states

  double Ec = shift;
  double Ev = shift;

  double new_shift = shift;

  bool foundall = true;

  // did we find all electron eigenstates?
  int n_eig = _upt_solver_options.n_vb;
  for ( ; n_eig < n_states; n_eig++)
  {
    if (_solution[n_eig].eigen_vector.size() == 0)
    {
      foundall = false;
      break;
    }
    else
      new_shift = _solution[n_eig].eigen_energy - 1e-5;
  }

  if (n_eig == _upt_solver_options.n_vb)
  {
    // in this case we found no electron state at all, so set shift near Ec
    new_shift = _upt_solver_options.guess_cb;
  }
  // if not all are found, look for so many electron states:
  number_of_eigenstates = _upt_solver_options.n_cb - (n_eig - _upt_solver_options.n_vb) + 1;

  if (foundall)
  {
    new_shift = shift;

    // did we find all hole eigenstates?
    for (n_eig = static_cast<int>(_upt_solver_options.n_vb) - 1; n_eig >= 0; n_eig--)
    {
      if (_solution[n_eig].eigen_vector.size() == 0)
      {
        foundall = false;
        break;
      }
      else
        new_shift = _solution[n_eig].eigen_energy + 1e-5;
    }

    if (n_eig == static_cast<int>(_upt_solver_options.n_vb) - 1)
    {
      // in this case we found no state at all, so set shift near Ev
      new_shift = _upt_solver_options.guess_vb;

    }
    number_of_eigenstates = n_eig + 1;

  }


  return(make_pair(number_of_eigenstates, new_shift));

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
      std::cout<< "("+get_name()+") including band shifts" <<std::endl;
      add_band_shifts();
    }

    if(_upt_options.potential_flag)
    {
      std::cout<< "("+get_name()+") passing potential" <<std::endl;
      add_pot_shifts();
    }

    write_shifts();
    
    if(_upt_options.d_states_correction)
    {
      std::cout<< "("+get_name()+") passing strain" <<std::endl;
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
  Messages::info("("+get_name()+") compute atom-projected charges"); 

  _el_atomic_charges.resize(_N_without_H, 0.0);
  _hl_atomic_charges.resize(_N_without_H, 0.0);

  compute_atomic_charges("el", _el_atomic_charges);
  compute_atomic_charges("hl", _hl_atomic_charges);

  //Print for debug charges on atoms
  Messages::info("("+get_name()+") print atom-projected charges on files"); 

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
    Messages::info("("+get_name()+") write wave functions on files");
    inst->write_states();
  }

  if (_upt_options.pdos_flag == true)
  {
  compute_pdos();
  }

}
	


//-------------------------------------------------------------------------

void ETB::parse_options(void)
{
  //std::cout << "("+get_name()+") parse_options() begin...";

  const ModelOptions& solopts = get_solver_options();
  

  // get kpoint
  libMesh::RealVectorValue k_vec;
  get_option("k_vector", k_vec);
  set_k_point(k_vec);

  // get kpoint as parameter
  Point kp(get_k_point(true));
  // the following seems problematic, since the upt k-point may become inconsistent with that
  // in the base class?
  //get_parameter("k_x", _upt_options.k_point[0], kp(0) );
  //get_parameter("k_y", _upt_options.k_point[1], kp(1) );
  //get_parameter("k_z", _upt_options.k_point[2], kp(2) );
  _upt_options.k_point[0] = kp(0);
  _upt_options.k_point[1] = kp(1);
  _upt_options.k_point[2] = kp(2);


  _upg_filename = get_option("upg_filename", "none");

  _upt_options.verbose = get_option("verbose", SimulationOptions::verbose());

  _upt_options.etb_dataset = get_option("dataset","");
  _upt_options.max_TB_order = get_option("max_TB_order", 2);

  std::string sparse_fmt = get_option("sparse_format", "full");
  sparse_fmt.copy(_upt_options.sparse_fmt, sparse_fmt.size() );

  _upt_options.check_bondmap = get_option("check_bondmap", false);

  _upt_options.relat_flag = get_option("relativistic", true);

  _upt_options.pdos_flag = get_option("write_pdos", false);

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
  string passivation_model = get_option("passivation_model","hybrid");
  if ( passivation_model == "hybrid" )
      {_upt_options.hybrid_passivation = true;}

  //---------------------------------------------------------------------------------------
 
  _upt_options.assemble_H = get_option("assemble_hamiltonian",true);


 
  // Solver options: "upt_lanczos"  
  _upt_solver_options.solver = solopts.get_option("solver", "lanczos");
  string solver_type = solopts.get_option("solver_type", "cpu");
  if ( solver_type == "cpu") _upt_solver_options.solver_flag = 0;
  if ( solver_type == "gpu") _upt_solver_options.solver_flag = 1;
  if ( solver_type == "gpu-split") _upt_solver_options.solver_flag = 2;

  _upt_solver_options.n_vb = get_option("num_valence_eigenvalues", 0);
  _upt_solver_options.n_vb =
      get_option("num_hole_states", _upt_solver_options.n_vb);
  _upt_solver_options.n_vb =
      solopts.get_option("num_valence_eigenvalues", _upt_solver_options.n_vb);
  _upt_solver_options.n_vb =
      solopts.get_option("num_hole_states", _upt_solver_options.n_vb);
  
  //if (_upt_solver_options.n_vb%2==1) _upt_solver_options.n_vb += 1;
  
  _upt_solver_options.n_cb =  get_option("num_conduction_eigenvalues", 0);
  _upt_solver_options.n_cb =
      get_option("num_electron_states", _upt_solver_options.n_cb);
  _upt_solver_options.n_cb =
      solopts.get_option("num_conduction_eigenvalues", _upt_solver_options.n_cb);
  _upt_solver_options.n_cb =
      solopts.get_option("num_electron_states", _upt_solver_options.n_cb);

  //if (_upt_solver_options.n_cb%2==1) _upt_solver_options.n_cb += 1;
  
  _upt_solver_options.min_iter =  solopts.get_option("min_iter", 30);
  _upt_solver_options.long_iter =  solopts.get_option("long_iter", 32);
  _upt_solver_options.max_iter =  solopts.get_option("max_iter", 10000);

  bool flag = solopts.get_option("dynamic_search",true);
  if (flag) _upt_solver_options.dynamic = 1;
  _upt_solver_options.bitoff = solopts.get_option("dynamic_offset", 0.1);

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

  // da togliere e leggere dal database: shift della banda di valenza (che e` 0)
  //_upt_options.vb_shift = options.get_option("vb_shift", 0.0);

  _upt_solver_options.fast_tol =  solopts.get_option("fast_tolerance", 1e-1);
  _upt_solver_options.long_tol =  solopts.get_option("long_tolerance", 1e-9);
  _upt_solver_options.ort_tol =  solopts.get_option("orthogonality_tolerance", 1e-6);

  //Get projection_length for quantum charge projection (Ang)
  _upt_options.projection_length = get_option("projection_length", 5.0);

  _upt_solver_options.guess_vb =
      solopts.get_option("guess_valence", _upt_solver_options.guess_vb);
  _upt_solver_options.guess_cb =
      solopts.get_option("guess_conduction", _upt_solver_options.guess_cb);

   //parse these to avoid warnings
   get_option("states_file", "none");
   get_option("load_file", "none");
   get_option("load_path", "none");

  //std::cout << "Projection lenght set to " <<  _upt_options.projection_length << std::endl;
  //std::cout << "done" << std::endl;

  if (plot_solution("eDensity") ||
      plot_solution("hDensity"))
    _upt_options.compute_densities = true;
}

//-------------------------------------------------------------------------

void
ETB::print_upt_options(void)
{

  std::cout << "("+get_name()+") UPTIGHT_OPTIONS: " << std::endl;

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
        if (atom[j].get_elem() != NULL)
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
      if (atom[i].get_elem() == NULL)
        _band_shift[i] = _band_shift[b_map[i][0]];
  }
  else
  {
    for (unsigned int i = 0; i < N; i++)
    {

      if (atom[i].get_elem() == NULL)
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

      // 6/5/2016 I believe one should take strain in crystal coordinates, here
      _strain_int.get_crystal_strain(el, el->centroid(), epsilon);

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

  for(j=0; j<N_atoms_wo_H; j++)
    qmat[j] = 0.0;

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

	// TODO it may be better to change sign of energies for holes, instead of doing the
	// 1-pop thing afterwards
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
ETB::do_project_to_primitive_cell(const vector<eigen_problem_solution>& a,
    const vector<Point>& kpoints,
    vector<vector<Complex>>& weights) const
{

  int n_k = kpoints.size();

  int n_states = a.size();
  vector<Complex> value(n_states, 0.0);


  const std::vector<Atom>& atom = get_atomistic_structure()->get_structure_atoms();
  size_t N = get_atomistic_structure()->get_N_without_H();

  int n_spin = _upt_options.relat_flag ? 2 : 1;

  vector<int> orbitals_i;
  vector<int> orbitals_j;

  size_t id_i = 0;

  for (size_t i = 0; i < N; i++)
  {
    // WARNING: _ion_num_orbitals considers the two spins!
    orbitals_i.resize(_ion_num_orbitals[i] / n_spin);
    inst->get_ion_orbitals(i+1, orbitals_i);

    size_t id_j = id_i;
    for (size_t j = i; j < N; j++)
    {
      if (atom[i].get_label() == atom[j].get_label())
      {

        orbitals_j.resize(_ion_num_orbitals[j] / n_spin);
        inst->get_ion_orbitals(j+1, orbitals_j);

        for (auto&& v : value)
          v = 0.0;

        for (size_t orb_i = 0; orb_i < orbitals_i.size(); ++orb_i)
        {
          for (size_t orb_j = 0; orb_j < orbitals_j.size(); ++orb_j)
          {
            if (orbitals_i[orb_i] == orbitals_j[orb_j])
            {
              for (int s = 0; s < n_states; ++s)
              {
                value[s] += a[s].eigen_vector[id_i + orb_i] *
                    libmesh_conj(a[s].eigen_vector[id_j + orb_j]);

                if (_upt_options.relat_flag)
                  value[s] += a[s].eigen_vector[id_i + orb_i + orbitals_i.size()] *
                  libmesh_conj(a[s].eigen_vector[id_j + orb_j + orbitals_j.size()]);
              }
            }
          }
        }

        Point d(atom[i].get_position());
        d.subtract(atom[j].get_position());

        double mul = (j == i) ? 1 : 2;

        for (int k = 0; k < n_k; ++k)
        {
          // 0.1 is for conversion from Angstrom to nm
          Complex arg = 0.1*Complex(0.0, kpoints[k]*d);
          Complex phase = std::exp(-arg);

          for (int s = 0; s < n_states; ++s)
            weights[k][s] += mul*libmesh_real(value[s] * phase);
        }
      }

      id_j += _ion_num_orbitals[j];
    }

    id_i += _ion_num_orbitals[i];
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

void ETB::compute_pdos(void)
{
  string out_path = get_output_directory();
  const std::vector<Atom>& atom = get_atomistic_structure()->get_structure_atoms();
  unsigned int n = _solution.size();
  unsigned int N_atoms_wo_H = get_atomistic_structure()->get_N_without_H();

  unsigned int k = 0;
  unsigned int k_at = 0;

  double s = 0.01;
  unsigned int step = 1000;
  std::vector<double> atom_pdos(step+1);

  double emin = _solution[0].eigen_energy;
  double emax = _solution[0].eigen_energy;
  for (unsigned int i = 1; i < n; i++)
  {
    if(_solution[i].eigen_energy < emin){emin = _solution[i].eigen_energy;}
    if(_solution[i].eigen_energy > emax){emax = _solution[i].eigen_energy;}
  }
  double de = (emax-emin)/step;

  for (unsigned int j = 0; j < N_atoms_wo_H; j++)
  {

    ofstream tb_pdos(out_path + "/pdos" + to_string(j+1) + ".dat");

    tb_pdos << "#eV   pdos   " << atom[j].get_specie() << endl;

    for(unsigned int i = 0; i <= step; i++){atom_pdos[i] = 0.0; }

    for (unsigned int i = 0; i < n; i++)
    {

      double atom_sum = 0.0;

      for (k = k_at; k < k_at + _ion_num_orbitals[j]; k++)
      {
        atom_sum += std::norm(_solution[i].eigen_vector[k]);
      }

      for (unsigned int l = 0; l <= step; l++)
      {
        atom_pdos[l] += atom_sum*1/(sqrt(2*3.141592653589793*s*s))*exp(-0.5*(((emin + l*de)-_solution[i].eigen_energy)/s)*(((emin + l*de)-_solution[i].eigen_energy)/s));
      }

    }

    for (unsigned int l = 0; l <= step; l++)
    {
      tb_pdos << emin + l*de << "  " << atom_pdos[l] << endl;
    }

    k_at = k;

  }
}

void ETB::get_valence_band_shifts(const Material* mat)
{
  if (mat->is_alloy())
  {
    const Alloy* alloy = dynamic_cast<const Alloy*>(mat);
    vector<const Material*> mats(2);
    mats[0] = alloy->get_component_A();
    mats[1] = alloy->get_component_B();
    for (int i = 0; i < mats.size(); ++i)
      get_valence_band_shifts(mats[i]);
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
          get_valence_band_shifts(mat);
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

    Messages::info("("+get_name()+") get potential from "+_upt_options.potential_sim);

    string band_src("");
    if (edge=="Ec")
      band_src = get_option("cb_edge", _dd_int->get_name() + ".Ec");
    else
      band_src = get_option("vb_edge", _dd_int->get_name() + ".Ev");

    SolutionProvider cb = find_solution_provider(band_src);

    if ((cb.second == INVALID_ID))
      throw RuntimeException("Simulation \'" + _dd_int->get_name() +
          "\' lacks band edge solution variables.");

    vector<double> edges(8,0.0);

    MeshBase::const_element_iterator it = this->active_local_elements_begin();
    const MeshBase::const_element_iterator end = this->active_local_elements_end();

    for ( ; it != end; ++it)
    { 
      const Elem* elem = *it;
      vector<Point> p(elem->n_nodes());

      for (size_t i = 0; i < elem->n_nodes(); ++i)
        p[i] = elem->point(i);

      cb.first->get_solution(elem, cb.second, edges, p);

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


  if (edge == "Ec")
    this->get_communicator().min(band_edge);
  else if (edge == "Ev")
    this->get_communicator().max(band_edge);

  return band_edge;

}


void
ETB::get_band_extrema(double& cb_min, double& vb_max)
{
  if (_upt_options.potential_flag)
  {
    Messages::info("("+get_name()+") get band edges from model: "+_upt_options.potential_sim);

    string cb_src = get_option("cb_edge", _dd_int->get_name() + ".Ec");
    string vb_src = get_option("vb_edge", _dd_int->get_name() + ".Ev");

    SolutionProvider cb = find_solution_provider(cb_src);
    SolutionProvider vb = find_solution_provider(vb_src);

    if ((cb.second == INVALID_ID) || (vb.second == INVALID_ID))
      throw RuntimeException("Simulation \'" + _dd_int->get_name() +
          "\' lacks band edge solution variables.");

    vb_max = -numeric_limits<double>::max();
    cb_min = numeric_limits<double>::max();

    //map<ID, vector<double> > bandedges;
    vector<double> cb_edges(8,0.0);
    vector<double> vb_edges(8,0.0);

    MeshBase::const_element_iterator it = this->active_local_elements_begin();
    const MeshBase::const_element_iterator end = this->active_local_elements_end();

    for ( ; it != end; ++it)
    { 
      const Elem* elem = *it;
      vector<Point> p(elem->n_nodes());

      for (size_t i = 0; i < elem->n_nodes(); ++i)
        p[i] = elem->point(i);

      cb.first->get_solution(elem, cb.second, cb_edges, p);
      vb.first->get_solution(elem, vb.second, vb_edges, p);

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

  this->get_communicator().min(cb_min);
  this->get_communicator().max(vb_max);

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

  bool do_edens = values.count(eDensity);
  bool do_hdens = values.count(hDensity);

  if (do_edens || do_hdens)
  {
    TiberLinearSystem& qdens_sys = get_equation_system<TiberLinearSystem>(0);
    NumericVector<Number>& qdens = *qdens_sys.solution;
    unsigned int dim = get_mesh().mesh_dimension();

    FEType fe_type = qdens_sys.variable_type(0);
    UniquePtr<FEBase> fe(build_finite_element(dim, fe_type));
    const vector<vector<Real> >& phi = fe->get_phi();

    fe->reinit(elem, &p);

    DofMap& dof_map = qdens_sys.get_dof_map();
    std::vector<unsigned int> dof_indices_el, dof_indices_hl;
    dof_map.dof_indices(elem, dof_indices_el, 0);
    dof_map.dof_indices(elem, dof_indices_hl, 1);
    unsigned int n_dofs = phi.size();

    for (unsigned int n = 0; n < np; n++)
    {
      double value_e = 0;
      double value_h = 0;

      // the phi^2 factor comes from the fact that the more correct interpolation is the square
      // of the basis function, because the probability densities are the square of the states
      // NOTE: maybe one should check if this gives really a better result
      // NOTE: 2011-12-01 the above turned out to be wrong: phi is used only as linear interp. !
      for (unsigned int i = 0; i < n_dofs; i++)
      {
         value_e += phi[i][n] * qdens(dof_indices_el[i]);
         value_h += phi[i][n] * qdens(dof_indices_hl[i]);
      }

      if (do_edens)
        values[eDensity][n] = value_e;
      if (do_hdens)
        values[hDensity][n] = value_h;
    }

  }

  /*
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
  */

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

  //Get point coordinate as physical_point
  Point phys_p;

  /*
  if (values.count(ElQuantumDensityNodes))
    {
    for (unsigned int n = 0; n < np; n++)
      {
      if (_dim == 3)
        {
        phys_p = libMesh::FE< 3, libMeshEnums::LAGRANGE>::map(elem, p[n]);
        values[ElQuantumDensityNodes][n] = build_rho3d(_el_atomic_charges, elem, phys_p);
        }
      else if (_dim == 2)
        {
        phys_p = libMesh::FE< 2, libMeshEnums::LAGRANGE>::map(elem, p[n]);
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
        phys_p = libMesh::FE< 3, libMeshEnums::LAGRANGE>::map(elem, p[n]);
        values[HlQuantumDensityNodes][n] = build_rho3d(_hl_atomic_charges, elem, phys_p);
        }
      else if (_dim == 2)
        {
        phys_p = libMesh::FE< 2, libMeshEnums::LAGRANGE>::map(elem, p[n]);
        values[HlQuantumDensityNodes][n] = build_rho2d(_hl_atomic_charges, elem, phys_p);
        }
      else if (_dim == 1)
        {
        throw ModelErrorException("HlQuantumNodes is not supported in 1D "
            "calculations. Use HlQuantumDensity (CELL variable)");
        }

      }
    }
  */

  if (values.count(MeshStatesNodes))
    {
    for (unsigned int n = 0; n < np; n++)
      {

      for (unsigned int i = 0; i < _solution_size; i++)
        {
        if (_dim == 3)
          {
          phys_p = libMesh::FE< 3, libMeshEnums::LAGRANGE>::map(elem, p[n]);
          values[MeshStatesNodes][_solution_size * n + i] = 
            build_rho3d(_eigenvector_mag[i], elem, phys_p);
          }
        else if (_dim == 2)
          {
          phys_p = libMesh::FE< 2, libMeshEnums::LAGRANGE>::map(elem, p[n]);
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
  //std::cout<<"("+get_name()+") get neigh_atoms "<<endl;
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

  //declare_solution(ElQuantumDensity, REAL, CELL, "q"+units);
  //declare_solution(HlQuantumDensity, REAL, CELL, "q"+units);
  declare_solution(eDensity, REAL, NODES, "q"+units);
  declare_solution(hDensity, REAL, NODES, "q"+units);
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

  // get the global matrix size
  int global_size = size_matrix;
  this->get_solver_communicator().sum(global_size);
  
  vector<int> non_zeros_number(size_matrix);
  vector<int> offdiag_nnz(0);
  int nnz = 0;

  for (int row = 0 ; row < size_matrix; row++)	
  {
    // nrow + 1 to get Fortran indexing
    non_zeros_number[row] = inst->get_H_row_size(row + 1);
    nnz += non_zeros_number[row];
  }


  EigenSolver::preallocate_matrix('H', global_size, size_matrix, non_zeros_number, offdiag_nnz);


  //----------------------------------------------------------------------------------------------------//
 
  //write data of columns in each row
  
  for (int row = 0 ; row < size_matrix; row++)
  {

      // nrow + 1 to get Fortran indexing
      int n_cols = inst->get_H_row_size(row + 1);

      vector<unsigned int> column_vector(n_cols);
      vector<Complex> row_values(n_cols);
 
      inst->get_H_row(row + 1, reinterpret_cast<int*>(column_vector.data()), row_values.data());

      EigenSolver::insert_matrix_row('H', row, column_vector, row_values, 1);
  }

  EigenSolver::finalize_matrix_assembly('H');
  
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
ETB::get_H_csr(std::vector<Complex>& A,
               std::vector<int>& JA,
               std::vector<int>& IA,
               const std::vector<unsigned int>& perm) const
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
    this->set_atomistic_structure(get_environment().get_device().get_atomistic_structure(name));
    if (get_atomistic_structure() == NULL)
      throw InitFailedException("No atomistic structure \'" + name + "\' found "
          "for simulation \'" + get_name());

    if (!includes_regions(get_atomistic_structure()->get_IDset()))
    { 
       Messages::warning("Module will restrict the atomistic structure '"+name+"'");      
       Messages::info("Note: the module will work on a copy of the original structure");      


       // TODO Warning: this looks like producing potential memory leakage
       AtomisticStructure* new_as = AtomisticStructure::create(*get_atomistic_structure());

       set_atomistic_structure(new_as);
    
    }
  }
}


pair<double, double>
ETB::project_densities(const Elem* elem, const Point& point, double cutoff)
{
  pair<double, double> densities(0, 0);

  AtomisticStructure* as = get_atomistic_structure();


  if (get_option("new_projection", false))
  {
    // TODO WARNING assuming cutoff in A, but find_nearest_atom expects nm
    int index = as->find_nearest_atom(elem, point, 2 * cutoff / 10.0);

    if (index >= 0)
    {

      const BondMap& bondmap = as->get_bond_map();
      int coordination = bondmap[index].size();

      if (coordination == 4)
      {
        BondMap::Translation tr = bondmap.get_translation();
        // take the first as origin
        Point o(as->get_structure_atom(bondmap[index][0]).get_position()
            + tr[index][0]);
        Point a(as->get_structure_atom(bondmap[index][1]).get_position()
            + tr[index][1] - o);
        Point b(as->get_structure_atom(bondmap[index][2]).get_position()
            + tr[index][2] - o);
        Point c(as->get_structure_atom(bondmap[index][3]).get_position()
            + tr[index][3] - o);

        // volume of the tetrahedron
        double vol = 1.0/6.0 * a * b.cross(c);
        double normalization = 0.5* 1e24 / fabs(vol);
        //double normalization = 1e24;

        densities.first += _el_atomic_charges[index] * normalization;
        densities.second += _hl_atomic_charges[index] * normalization;

        for (unsigned int k = 0; k < 4; ++k)
        {
          unsigned int atom_id = bondmap[index][k];
          int coordination = bondmap[atom_id].size();
          if (coordination == 4)
          {
            densities.first += _el_atomic_charges[atom_id] * normalization / coordination;
            densities.second += _hl_atomic_charges[atom_id] * normalization / coordination;
          }
        }
      }
    }
  }
  else
  {

    const double scale = get_atomistic_structure()->get_scale();
    const double sigma = cutoff;
    const double sigma2 = 2.0*sigma*sigma;

    // the point in Angstrom
    Point coord(point);
    coord *= get_mesh_units() / 1e-10;

    // normalization/scale factor
    double normalization = 1.0;
    switch (_dim)
    {
//      case 1:
//        normalization = 1e8;
//        break;

//      case 2:
//        normalization = 1e16 / (2.0 * M_PI * sigma);
//        break;

      default:
      {
        double tmp = 1 / (2.0 * M_PI * sigma * sigma);
        normalization = 1e24 * sqrt(tmp * tmp * tmp);
        break;
      }
    }
    //normalization = 3.0 / (4 * M_PI * sigma * sigma * sigma) * 1e24;

    // this approach does not need the elem_atoms map, but then it needs to be able
    // to find atoms not only in spheres.
    ///*
    int index = as->find_nearest_atom(elem, point, 2 * cutoff / 10.0);
    AtomisticBasis::neighbor_iterator it(
        get_atomistic_structure()->neighbors_begin(index, cutoff));
    const AtomisticBasis::neighbor_iterator end(
        get_atomistic_structure()->neighbors_end(index, cutoff));

    // the charge of the nearest atom
    //densities.first += _el_atomic_charges[index] * normalization;
    //densities.second += _hl_atomic_charges[index] * normalization;

    for (; it != end; ++it)
    {
      const Atom* atom = *it;
      unsigned int atom_id = it.atom_index();

      //cerr << " " << atom_id << endl;
      Point atom_pos(atom->get_position() + it.atom_translation());
    //*/

    /*
    const std::vector<unsigned int>& atoms = get_elem_atoms(elem->id());
    for (int i = 0; i < atoms.size(); ++i)
    {
      unsigned int atom_id = atoms[i];
      const Atom* atom = &get_atomistic_structure()->get_structure_atom(atom_id);
      Point atom_pos(atom->get_position());
    */

      Point delta_r = coord - atom_pos;

      switch (_dim)
      {
//        case 1:
          // set dy = 0, the rest is the same as in 3D
//          delta_r(1) = 0.0;

//        case 2:
          // set dz = 0, the rest is the same as in 3D
//          delta_r(2) = 0.0;

        default:
        {
          double factor = normalization * exp(-delta_r.norm_sq() / sigma2);
          densities.first += _el_atomic_charges[atom_id] * factor;
          densities.second += _hl_atomic_charges[atom_id] * factor;
          //densities.first += _el_atomic_charges[atom_id] * normalization / coordination;
          //densities.second += _hl_atomic_charges[atom_id] * normalization / coordination;
          break;
        }
      }
    }
   //cerr << densities.first << endl;

  }

  return densities;
}
