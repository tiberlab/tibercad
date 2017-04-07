#include "KmcTUM.h"
#include "TiberCad.h"
#include "TiberModule.h"
#include "InitFailedException.h"
#include "MeshUtils.h"
#include "Constants.h"
#include "mesh_generation.h" 
#include "Database.h"
#include "SimulationEnvironment.h"
#include "SimulationOptions.h"
#include "Messages.h"
#include "node.h"
#include <sstream> 
#include <fstream> 
#include <cstdlib>
#include "gmsh_io.h"
#include <boost/filesystem.hpp>
#include "boost/process.hpp"

using namespace std;
namespace bp = boost::process;
namespace bpi = boost::process::initializers;
namespace fs = boost::filesystem;

KmcTUM::KmcTUM(const ModelOptions& options) :
  KmcInterface(options)
{
   has_solution_vector(false);
}

KmcTUM*
KmcTUM::create(const ModelOptions& opt)
{
  return new KmcTUM(opt);
}


void
KmcTUM::do_setup_solution_variables(void)
{
  declare_solution(elDensity, REAL, CELL, "1/cm^3");
  declare_solution(hlDensity, REAL, CELL, "1/cm^3");
  //declare_solution(xDensity, REAL, CELL, "1/cm^3");
  declare_solution(eCurrentDensity, REAL, CELL, "A/cm^2");
  declare_solution(hCurrentDensity, REAL, CELL, "A/cm^2");
  //declare_solution(Recombination, REAL, CELL, "s^-1 cm^-3");
  //declare_solution(Generation, REAL, CELL, "s^-1 cm^-3");
}



void
KmcTUM::do_init()
{ 
  KmcInterface::do_init();
  
  parse_options();  
 
  _eldensity.resize(grid.num_elements()); 
  _hldensity.resize(grid.num_elements()); 
  _potential.resize(grid.num_elements());

  
  // get the poisson bounding box in order to define the buffer distances  
  SimulationInterface* poisson_sim = _pot_sol.first;

  pair<Point, Point> bbox = poisson_sim->get_environment().get_bounding_box(true);
 
  double dz = grid.grid_step(0);
  unsigned int nz = grid.num_elements(2);

  _buffer_bottom = _p0(2) - bbox.first(2) + dz/2.0;
  _buffer_top = bbox.second(2) - _p0(2) - nz*dz + dz/2.0;
  
  cout<<"BUFFERS: "<<_buffer_bottom<<" "<<_buffer_top<<endl;

  //cout<<"(kmc) eq sys size: "<<grid.num_elements()<<endl;
  // setup a linear system for interpolation purposes
  //ID id = create_equation_system("linear","");
  //_sys = &get_equation_system<TiberLinearSystem>(id);
  //cout<<"(kmc) eq sys initialized"<<endl;

}
 
void
KmcTUM::setup_mesh(void)
{
  // setup standard mesh
  SimulationInterface::setup_mesh();

  // 
  // The grid points of the kmc code, 'o' are the centroids
  // of the mesh elements of this module (+ are the nodes)
  //
  //           XY                        Z        
  //    0 1           N-1        N   o o o o o          
  //   +-+-+-+-+-+-+-+-+            +-+-+-+-+-+   
  //   |o|o|o|o|o|o|o|o|o       N-1 |o|o|o|o|o|   
  //   +-+-+-+-+-+-+-+-+            +-+-+-+-+-+ 
  //   |o|o|o|o|o|o|o|o|o           |o|o|o|o|o| 
  //   +-+-+-+-+-+-+-+-+            +-+-+-+-+-+ 
  //   |o|o|o|o|o|o|o|o|o           |o|o|o|o|o| 
  //   +-+-+-+-+-+-+-+-+            +-+-+-+-+-+ 
  //   |o|o|o|o|o|o|o|o|o        1  |o|o|o|o|o|
  //   +-+-+-+-+-+-+-+-+            +-+-+-+-+-+
  //   0 1             N         0   o o o o o   
  //                           
 
  // we get the bounding box of the kmc region
  pair<Point, Point> bbox = get_environment().get_bounding_box(true);

  double dl[3]; 
  unsigned int nx[3];
  // read-in grid parameters for kmc
  double gs = get_option("grid_step", 1.00); // in nanometers
  dl[0] = get_option("grid_step_a",gs); // in nanometers
  dl[1] = get_option("grid_step_b",dl[0]); // in nanometers
  dl[2] = get_option("grid_step_c",dl[1]); // in nanometers
  
  Point p1(bbox.first);
  Point p2(bbox.second);
  Point sides=p2-p1;
 
  nx[0] = static_cast<int>(sides(0)/dl[0]); 
  nx[1] = static_cast<int>(sides(1)/dl[1]); 
  nx[2] = static_cast<int>(sides(2)/dl[2]) + 1;  
 
  Point dr(dl[0]/2.0, dl[1]/2.0, -dl[2]/2.0); 
  _p0 = p1+dr; // store the origin shifted 

  //cout<<"region bbox: "<<p1<<" : "<<p2<<endl;
  //cout<<"region sides: "<<sides<<p2<<endl;
  //cout<<"grid setup: "<<nx[0]<<"x"<<nx[1]<<"x"<<nx[2]<<endl;

  grid.setup(p1,p2,nx[0],nx[1],nx[2]);

  /* DEFINITION OF MESH IS GOOD FOR INTERPOLATION
  MeshTools::Generation::build_cube(*box_mesh, 
			       nx[0],
			       nx[1],
			       nx[2],
		         p1(0), p2(0),
				     p1(1), p2(1),
				     p1(2), p2(2),
				     HEX8);

  _es=new(*box_mesh);
  System system = _es.add_system<ExplicitSystem>("boxmesh");
  system->add_variable("edens", FIRST);
  system->add_variable("hdens", FIRST);
  _es.init();
  */
}

void
KmcTUM::do_reinit()
{
  // Read Poisson potential and print on a file
  {
    SimulationInterface* poisson_sim = _pot_sol.first;
    ID sol_id = _pot_sol.second;
 
    MeshUtils::GridMapper& mapper =
        MeshUtils::GridMapper::get_mapper(poisson_sim->get_mesh());
 
    for (unsigned int i = 0; i < grid.num_elements(); ++i)
    {
      Point p(grid.get_centroid(i));
      double value;
      const Elem* elem = mapper.get_element(p);
      if ((elem != NULL) && poisson_sim->get_solution(elem, sol_id, value, p))
      {
        _potential[i] = value; // V
      }
    }
  }

  // Read chemical potentials at contacts
  {
    SimulationEnvironment* env = &get_environment();
    SimulationEnvironment::BoundaryIterator it = env->boundaries_begin();
    SimulationEnvironment::BoundaryIterator itend = env->boundaries_end();
    
    map<Boundary*, double> zb;

    unsigned int nx = grid.num_elements(0);
    unsigned int ny = grid.num_elements(1);
    unsigned int nz = grid.num_elements(2);
    double dx = grid.grid_step(0);
    double dy = grid.grid_step(0);
    double dz = grid.grid_step(0);
    Point pt;
    Point pb(0.0);
    double phi, efermi, hfermi;  
    double av_phi, av_efermi, av_hfermi;
    
    for (; it != itend; ++it)
    {
       SimulationEnvironment::BoundarySideIterator sit = 
                       env->boundary_sides_begin((*it)->get_name());

       const Node* node = (*sit).elem()->side((*sit).side())->node_ptr(0);
       zb[*it] = (*node)(2);

       pb = _p0;
       pb(2) = zb[*it]; 

       av_phi = 0.0; av_efermi= 0.0; av_hfermi=0.0;
      
       for (unsigned int j=0; j < ny; j++)
       {
         for (unsigned int i=0; i < nx; i++)
         {
            Point p(dx*i,dy*j,0);
            pt = pb + p;          
            get_boundary_potentials(*it, pt, phi, efermi, hfermi);
            av_phi += phi;
            av_efermi += efermi;
            av_hfermi += hfermi;
         } 
       }
       av_phi = av_phi / (nx*ny);  
       av_efermi = av_efermi / (nx*ny);  
       av_hfermi = av_hfermi / (nx*ny);
   
       _mue[*it] = av_efermi;
       _muh[*it] = av_hfermi;
    }

    // The following is punk...  
    it = env->boundaries_begin();
    itend = ++it;
    it = env->boundaries_begin();

    if (zb[*it] < zb[*itend])  
    {
       _mue_cathode = _mue[*it]; 
       _muh_cathode = _muh[*it]; 
       _mue_anode = _mue[*itend]; 
       _muh_anode = _muh[*itend]; 
    }else{   
       _mue_cathode = _mue[*itend]; 
       _muh_cathode = _muh[*itend]; 
       _mue_anode = _mue[*it]; 
       _muh_anode = _muh[*it]; 
    }

  }
  


}

void
KmcTUM::do_solve()
{
  write_input();

  write_potential(); 

  fs::path exe(_exe_name);  
  fs::path lib("../lib");
  set<string> ld_lib_path;
  ld_lib_path.insert("LD_LIBRARY_PATH=" + (exe.branch_path()/lib).string() );
  vector<string> args;
  args.push_back(fs::absolute(exe).string());
  args.push_back(_kmc_scratch+"/input.txt");

  cout<<"Starting child process ... "<<endl;
 
  bp::child c = bp::execute(
      bpi::set_args(args),
      bpi::set_env(ld_lib_path),
      bpi::start_in_dir(".")
  );
  wait_for_exit(c);
  
  if (!read_density()) 
    Messages::error("density file error");  

}

void
KmcTUM::parse_options(void)
{

  ModelOptions &opts = get_solver_options();

  if (!opts.find_option("executable_name") )
     throw InitFailedException("User must specify 'executable_name' with full path"); 
  else
    _exe_name = opts.get_option("executable_name","osc");

  fs::path exe(_exe_name);
  if (!exists(exe) || !is_regular_file(exe))
     throw InitFailedException("Invalid executable to 'executable_name'");

  _kmc_scratch = opts.get_option("input_folder",get_scratch_directory());

  fs::path scratch(_kmc_scratch);
  if (exists(scratch) && !is_directory(scratch))
  {
     throw InitFailedException("Cannot create scratch folder: it seems file already exists");
  }
  else
  {
     fs::create_directory(scratch); 
  }
  
  // get a potential simulation
  // Format: "potential_simulation=poisson.Potential"
  string poisson = get_option("potential_simulation", "");
  if (poisson.empty())
     throw InitFailedException("potential_simulation must be defined");
 
  _pot_sol = this->find_solution_provider(poisson,"Potential");
  if (_pot_sol.first==NULL || _pot_sol.second == INVALID_ID)
     throw InitFailedException(poisson + " does not exist or does not provide potential");
      
  // get a simulation for Fermi levels
  // Format: "qfermi_simulation=driftdiffusion"
  // By default we use the same as potential but need to change variables
  // first we isolate the first token corresponding to the simulation name
  vector<string> tokens;
  Utils::tokenize(poisson, tokens, ".");
  string qfermi = get_option("qfermi_simulation", tokens[0]);
  _mue_sol = this->find_solution_provider(qfermi,"eQFermi");
  if (_mue_sol.first==NULL || _mue_sol.second == INVALID_ID)
     throw InitFailedException(qfermi + " does not exist or does not provide eQFermi levels");

  _muh_sol = this->find_solution_provider(qfermi,"hQFermi");
  if (_muh_sol.first==NULL || _muh_sol.second == INVALID_ID)
     throw InitFailedException(qfermi + " does not exist or does not provide hQFermi levels");
      

  _Nstep        = opts.get_option("Nstep",1000);
  _IOstep       = opts.get_option("IOstep",_Nstep);
  eps_r         = opts.get_option("eps_r",1.0);
  _create_blend = opts.get_option("create_blend",false);
  _volume_ratio = opts.get_option("blend_volume_ratio",1.0);
  _generation_steps= opts.get_option("blend_generation_steps",5000);

  if (!opts.find_option("blend_acceptor")) 
      throw InitFailedException("option 'blend_acceptor' not found");
  if (!opts.find_option("blend_donor")) 
      throw InitFailedException("option 'blend_donor' not found");

  // getting database default values for HOMO/LUMO of acceptor
  string blend1 = opts.get_option("blend_acceptor","none");
  Database db1(blend1, blend1+".dat");
  db1.set_section("valenceband");

  _HOMO_acceptor= opts.get_option("HOMO_acceptor",db1.get("E_v",0.0,true));
  db1.set_section("conductionband");
  _LUMO_acceptor= opts.get_option("LUMO_acceptor",db1.get("E_c",0.0,true));

  // getting database default values for HOMO/LUMO of donor 
  string blend2 = opts.get_option("blend_donor","none");
  Database db2(blend2, blend2+".dat");  
  db2.set_section("valenceband");
  _HOMO_donor   = opts.get_option("HOMO_donor   ",db2.get("E_v",0.0,true));
  db2.set_section("conductionband");
  _LUMO_donor   = opts.get_option("LUMO_donor   ",db2.get("E_c",0.0,true));

  _sigma        = opts.get_option("disorder",0.050);
  _x_v0         = opts.get_option("exciton_attempt_rate",2.0e11);
  _x_sep        = opts.get_option("exciton_separation_rate",2.0e13);
  _x_rec        = opts.get_option("exciton_recombination_rate",2.0e9);
  _el_v0        = opts.get_option("el_attempt_rate",3.0e11);
  _hl_v0        = opts.get_option("hl_attempt_rate",6.0e11);
  _el_inj       = opts.get_option("el_injection_rate",3.0e12);
  _hl_inj       = opts.get_option("hl_injection_rate",1.0e10); 
  _alpha        = opts.get_option("alpha",2.0);
  _el_hl_R      = opts.get_option("el_hl_recombination_rate",5.0e5);
  _el_coll      = opts.get_option("el_collection_rate",1.0e10);
  _hl_coll      = opts.get_option("hl_collection_rate",1.0e10); 
  _coll_el_cathode = opts.get_option("max_el_collection_cathode",20000);
  _coll_hl_cathode = opts.get_option("max_hl_collection_cathode",20000);
  _coll_el_anode   = opts.get_option("max_el_collection_anode",20000);
  _coll_hl_anode   = opts.get_option("max_hl_collection_anode",20000);
  _kmc_total_time = opts.get_option("simulation_time",1e-3);
  _kmc_stat_time = opts.get_option("time_interval_for_averges",1e-5);
  _const_gen_rate = opts.get_option("constant_generation_rate",0.0);

  // faster (I think) to use array than vector;
  //vector<bool> periodic(3,false);
  //get_option("periodic",periodic);

  //for (int i=0; i<3; i++)
  //  _periodic[i] = periodic[i];
  
  // contact work-functions: in general we need to read from drift-diffusion module 
  // needs to make a boundary iterator
  SimulationEnvironment* env = &get_environment();
  SimulationEnvironment::BoundaryIterator it = env->boundaries_begin();
  SimulationEnvironment::BoundaryIterator itend = env->boundaries_end();
  unsigned int count = 0;

  for (; it != itend; ++it)
  {
     count++;
     const ModelOptions& opts = (*it)->get_options();
     _mue[*it] = opts.get_option("qfermi_e", 0.0);

     _muh[*it] = opts.get_option("qfermi_h", 0.0);
  }
  if (count > 2){
    throw InitFailedException("Only two boundaries must be defined");
  } 

}


void 
KmcTUM::get_solution_secure(const Elem* elem,
                           map<ID, vector<double> >& values,
                           const vector<Point>& p)
{
  
   if (values.count(elDensity))
   { 
     for (unsigned int n=0; n < p.size(); n++)
     { 
       unsigned int i = grid.find_element(p[n]);
       values[elDensity][n] = _eldensity[i];  // 1/cm^3
     }
   }
   
   if (values.count(hlDensity))
   {  
     for (int n=0; n< p.size(); n++)
     { 
       unsigned int i = grid.find_element(p[n]);
       values[hlDensity][n] = _hldensity[i];  // 1/cm^3
     }
   }
  
  

}

double
KmcTUM::linear_interpolation(const Point& p)
{
  /*
   // search the 3 closest grid points 
   // and perform a linear interpolation
   unsigned int i = grid.find_element(p);
   Point p0 = grid.get_centroid(i);
   unsigned int k,l,m; 
   grid.element_to_index(i,k,l,m);
   
   int sgnx = p(0)-p0(0) > 0 ? +1 : -1;
   Point p1(p0);
   p1(0) += sgnx*grid.grid_step(0);
   k = mod(k+sgnx,grid.num_elements(0));

   int sgny = p(1)-p0(1) > 0 ? +1 : -1;
   Point p2(p0);
   p2(1) += sgny*grid.grid_step(1);
   l = mod(l+sgny,grid.num_elements(1));
   
   int sgnz = p(2)-p0(2) > 0 ? +1 : -1;
   Point p3(p0);
   p3(2) += sgnz*grid.grid_step(2);
   m = mod(m+sgnz,grid.num_elements(2));
  */
}       
       
void
KmcTUM::plot_globaldata(void)
{

}

void
KmcTUM::write_boundary_potential(void)
{
  string potfile = _kmc_scratch + "/potential.txt";
  fstream ff(potfile.c_str(), std::fstream::out);
  // needs to make a boundary iterator
  SimulationEnvironment* env = &get_environment();
  SimulationEnvironment::BoundaryIterator it = env->boundaries_begin();
  SimulationEnvironment::BoundaryIterator itend = env->boundaries_end();

  if(it == itend){Messages::error("Problem: you must define internal physical regions as BC");}
  

  unsigned int nx = grid.num_elements(0);
  unsigned int ny = grid.num_elements(1);
  unsigned int nz = grid.num_elements(2);
  double dx = grid.grid_step(0);
  double dy = grid.grid_step(0);
  double dz = grid.grid_step(0);
  Point pt;
  double phi, efermi, hfermi;  

  for (; it != itend; ++it)
  {
     cout<<"getting pot "<<(*it)->get_name()<<endl;
     ff<<"boundary: "<<(*it)->get_name()<<endl;
     ff<<"phi(x,y) x-fast y-slow "<<endl;
     for (unsigned int j=0; j < ny; j++)
     {
       for (unsigned int i=0; i < nx; i++)
       {
          Point p(dx*i,dy*j,0);
          pt = _p0 + p;           
          get_boundary_potentials(*it, pt, phi, efermi, hfermi);
        ff << phi << "  ";
     } 
     ff << endl;
   }  

}
ff.close();

}


void
KmcTUM::write_input(void)
{
  string inputpath = _kmc_scratch; 
  string outputpath = TiberCad::get_output_dir();
  string input_file = inputpath + "/input.txt";
  fstream ff(input_file.c_str(), std::fstream::out);

  ff<<"// input directory" << endl;
  ff<<"./"+inputpath+"/"<<endl;

  ff<<"// output directory" << endl;
  ff<<"./"+outputpath+"/"<<endl;
 
  ff<<"// load morphology or generate new (1 = load, 0 = generate)"<<endl;
  ff<<(_create_blend ? 0:1)<<endl;
  
  ff<<"// number of nodes in x-direction"<<endl;
  ff<<grid.num_elements(0)<<endl;
  ff<<"// number of nodes in y-direction"<<endl;
  ff<<grid.num_elements(1)<<endl;
  ff<<"// number of nodes in z-direction"<<endl;
  ff<<grid.num_elements(2)-1 <<endl;

  ff<<"// morphology generation: #steps"<<endl;
  ff<<_generation_steps<<endl;

  ff<<"// buffer layer top (nm):"<<endl;
  ff<<_buffer_top<<endl;
  ff<<"// buffer layer bottom (nm):"<<endl;
  ff<<_buffer_bottom<<endl;

  ff<<"// lattice parameter (nm)"<<endl;
  ff<<grid.grid_step(0)<<endl;

  ff<<"// HOMO of acceptor material (eV) "<<endl;
  ff<<_HOMO_acceptor<<endl; //(PCBM -6.1 to -6.8)
  ff<<"// LUMO of acceptor material (eV) "<<endl;
  ff<<_LUMO_acceptor<<endl; //(PCBM -3.7 to -4.3)
  ff<<"// HOMO of donor material (eV)"<<endl;
  ff<<_HOMO_donor<<endl; // (rr-P3HT -5.0 to -4.8) (PCBCTBT -5.0 to -5.3)
  ff<<"// LUMO of donor material (eV) "<<endl;
  ff<<_LUMO_donor<<endl; //(-3.0 to -2.9) (PCBCTBT -3.6)


  ff<<"// quasi fermi level e bottom (eV)"<<endl;
  ff<<_mue_cathode<<endl; 
  ff<<"// quasi fermi level h bottom (eV)"<<endl;
  ff<<_muh_cathode<<endl; 
  ff<<"// quasi fermi level e top (eV)"<<endl;
  ff<<_mue_anode<<endl; 
  ff<<"// quasi fermi level h top (eV)"<<endl;
  ff<<_muh_anode<<endl; 



  ff<<"// relative permittivity "<<endl;
  ff<<eps_r<<endl;

  ff<<"// correlation (1 = on, 0 = off)"<<endl;
  ff<<0<<endl;
  ff<<"// correlation radius (nm)"<<endl;
  ff<<3.<<endl;

  ff<<"//sigma_disorder (eV)"<<endl;
  ff<<_sigma<<endl;

  ff<<"// volume ratio acceptor:donor (1:x)"<<endl;
  ff<<_volume_ratio<<endl;

  ff<<"// minimum voltage (V)"<<endl;
  ff<<_bias<<endl;
  ff<<"// maximum voltage (V)"<<endl;
  ff<<_bias<<endl;
  ff<<"// voltage steps (V)"<<endl;
  ff<<1.0<<endl;

  ff<<"// simulated time [termination condition] (s)"<<endl;
  ff<<_kmc_total_time<<endl;
  ff<<"// delta t for current statistics evaluation (s)"<<endl;
  ff<<_kmc_stat_time<<endl;

  ff<<"// max number of electrons to be collected at top [termination] (#)"<<endl;
  ff<<_coll_el_cathode<<endl;
  ff<<"// max number of electrons to be collected at bottom [termination] (#)"<<endl;
  ff<<_coll_el_anode<<endl;
  ff<<"// max number of holes to be collected at top [termination] (#)"<<endl;
  ff<<_coll_hl_cathode<<endl;
  ff<<"// max number of holes to be collected at bottom [termination] (#)"<<endl;
  ff<<_coll_hl_anode<<endl;

  ff<<"// temperature (K)"<<endl;
  ff<<SimulationOptions::temperature<<endl;

  ff<<"// variable range hopping x-direction (# nodes)"<<endl;
  ff<<1<<endl;
  ff<<"// variable range hopping y-direction (# nodes)"<<endl;
  ff<<1<<endl;
  ff<<"// variable range hopping z-direction (# nodes)"<<endl;
  ff<<1<<endl;

  ff<<"// FLAG: dark current on/off (1 = on, 0 = off)"<<endl;
  ff<<1<<endl;
  ff<<"// FLAG: illumination on/off (1 = on, 0 = off)"<<endl;
  ff<<1<<endl;
  ff<<"// FLAG: Coulomb interaction on/off (1 = on, 0 = off)"<<endl;
  ff<<1<<endl;

  ff<<"// FLAG: kmcout on/off (1 = on, 0 = off) "<<endl;
  ff<<0<<endl; // attention: very large files (up to TBytes !)
  ff<<"// FLAG: statcurrout on/off (1 = on, 0 = off)"<<endl;
  ff<<0<<endl;
  ff<<"// FLAG: ivout on/off (1 = on, 0 = off)"<<endl;
  ff<<0<<endl;

  ff<<"// RATES: constant exciton generation rate (s^-1 nm^-3) "<<endl; 
  ff<<_const_gen_rate<<endl;   //set to 0. if external generation rates should be loaded

  ff<<"// RATES: exciton hopping rate (s^-1)"<<endl;
  ff<<_x_v0<<endl;
  ff<<"// RATES: exciton separation rate (s^-1)"<<endl;
  ff<<_x_sep<<endl;
  ff<<"// RATES: exciton decay rate (s^-1)"<<endl;
  ff<<_x_rec<<endl;
  ff<<"// RATES: electrons attempt-to-hop rate (s^-1)"<<endl;
  ff<<_el_v0<<endl;
  ff<<"// RATES: holes attempt-to-hop rate (s^-1)"<<endl;
  ff<<_hl_v0<<endl;

  ff<<"// injection rate prefactor"<<endl;
  ff<<1.<<endl; // better use max injection rates below and leave this at 1.0
  ff<<"// max injection rate electrons at cathode (s^-1)"<<endl;
  ff<<_el_inj<<endl; // (bottom, acceptor/Al interface) 
  ff<<"// max injection rate electrons at anode (s^-1)"<<endl;
  ff<<_el_inj<<endl; // (top, acceptor/PEDOT:PSS/ITO interface) 
  ff<<"// max injection rate holes at cathode (s^-1)"<<endl;
  ff<<_hl_inj<<endl; // (bottom, acceptor/Al interface)   
  ff<<"// max injection rate holes at anode (s^-1)"<<endl;
  ff<<_hl_inj<<endl; // (top, donor/PEDOT:PSS/ITO interface) 

  ff<<"// localization constant Miller-Abrahams (m^-1)"<<endl;
  ff<<_alpha<<endl;

  ff<<"// RATES: electron-hole recombination rate (s^-1)"<<endl;
  ff<<_el_hl_R<<endl;
  ff<<"// RATES: electron collection rate cathode (s^-1)"<<endl;
  ff<<_el_coll<<endl;
  ff<<"// RATES: electron collection rate anode (s^-1)"<<endl;
  ff<<_el_coll<<endl; //[change for EBL, e.g. PEDOT:PSS] 
  ff<<"// RATES: hole collection rate cathode (s^-1)"<<endl;
  ff<<_hl_coll<<endl; // [change for HBL]
  ff<<"// RATES hole collection rate anode (s^-1)"<<endl;
  ff<<_hl_coll<<endl;
  ff.close();
}

  
bool 
KmcTUM::read_density(void)
{

  string outputpath = TiberCad::get_output_dir();
  string input_file = outputpath + "/kmc-density.txt";
  fstream ff(input_file.c_str(), std::fstream::in);

  if (ff.fail()) return false;

  string tmp;
  for(unsigned int i=0; i < 5; i++)
    ff >> tmp;

  unsigned int k, l, m;
  double edens, hdens, pot;
  
  while (ff >> k >> l >> m >> edens >> hdens >> pot)
  {
     unsigned int i = grid.index_to_element(k,l,m);
     _eldensity[i] = edens;
     _hldensity[i] = hdens;
  }

  return true;

}

bool
KmcTUM::read_current(void)
{
  string outputpath = TiberCad::get_output_dir();
  string input_file = outputpath + "/kmc-current.txt";
  fstream ff(input_file.c_str(), std::fstream::in);
 
  if (ff.fail()) return false;

  string tmp;
  getline(ff, tmp);
  getline(ff, tmp);
  ff>> _ecurr_bottom;
  getline(ff,tmp);
  ff>> _hcurr_bottom;
  getline(ff,tmp);
  ff>> _ecurr_top;
  getline(ff,tmp);
  ff>> _hcurr_top;

  return true;

}
 //--------------------------------------------------- 
 //grid definition along Z 
 //            Z                               
 //    N   o o o o o    
 //       +-+-+-+-+-+   
 //   N-1 |o|o|o|o|o|   
 //       +-+-+-+-+-+   
 //       |o|o|o|o|o|   
 //       +-+-+-+-+-+   
 //       |o|o|o|o|o|    
 //       +-+-+-+-+-+    
 //    1  |o|o|o|o|o|     
 //       +-+-+-+-+-+       
 //    0   o o o o o 
 //--------------------------------------------------- 
 //
void
KmcTUM::write_potential(void)
{
   
  string inputpath = _kmc_scratch; 
  string input_file = inputpath + "/tcad-pot.txt";
  fstream ff(input_file.c_str(), std::fstream::out);

  unsigned int k,l,m;
  for (unsigned int i=0; i<grid.num_elements(); i++)
  {
     grid.element_to_index(i,k,l,m);
     if (m != 0 && m != grid.num_elements(2))
        ff<< k <<" "<< l <<" "<< m <<" "<< -(_potential[i]) <<endl;
  }

}
