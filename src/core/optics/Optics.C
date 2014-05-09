// $Id$

#include "Optics.h"
#include "EigenvalueProblem.h"
#include "SimulationInterface.h"
#include "SimulationOptions.h"
#include "AtomisticStructure.h"
#include "Messages.h"
#include "DataOutput.h"
#include "Constants.h"

using namespace std;
using namespace Constants; 


Optics::~Optics()
{
  delete(_energy_mesh);
}



Optics::Optics(const ModelOptions& options)
 : SimulationInterface(options)
{
  _initial_state_model = NULL;
  _final_state_model = NULL;
  _energy_mesh = NULL;
  _k_vector[0]=0.0;   _k_vector[1]=0.0;   _k_vector[2]=0.0; 
  _k_integration = NULL;

  has_solution_vector(false);
}

//==============================================//




void
Optics::do_setup_solution_variables(void)
{
  int dim = get_mesh().mesh_dimension();
  string unit("");
  if (dim == 0)
    unit = "/cm^3";
  if (dim == 1)
    unit = "/cm^2";
  else if (dim == 2)
    unit = "/cm";

  declare_solution(OpticalPower, REAL, GLOBAL, "W" + unit);
  declare_solution(Recombination, REAL, GLOBAL, "1/s" + unit);
}


void
Optics::get_solution_secure(map<ID, vector<double> >& values)
{

  map<ID, vector<double> >::iterator mapit(values.begin());
  const map<ID, vector<double> >::iterator mapend(values.end());
  for ( ; mapit != mapend; ++mapit)
  {
    ID id = mapit->first;

    switch (id)
    {
      case OpticalPower:
        values[id] = vector<double>(1, _total_power);
        break;

      case Recombination:
        values[id] = vector<double>(1, _recombination);
        break;

    }
  }
}



void Optics::parse_options()
{
   
  _initial_state_particle = "el";
  _final_state_particle = "hl";


  //k-vector
  RealVectorValue k_vec(0.0);
  get_parameter("k_vector",k_vec);
  for (short i = 0; i < 3; i++) _k_vector[i] = k_vec(i);
  

  std::string  job_name = get_option("job","matrix_elements");

  if (job_name == "matrix_elements")
    job = MATREL;
  else if (job_name == "bulk")
    job = BULKMATREL;
  else
    throw InitFailedException( "OpticsKP: Incorrect job: " + job_name);


  if (!plot_solution("optical_spectrum") && !plot_solution("optical_spectrum_k_0"))
      throw InitFailedException("(optics) requires a correct plot option");    

  if (plot_solution("optical_spectrum"))
    if (!get_options().has_submodel("k_integration")) 
      throw InitFailedException("(optics) requires a k-integration");

  
  if (plot_solution("optical_spectrum") && plot_solution("optical_spectrum_k_0") )
    throw InitFailedException("(optics) optical_spectrum incompatible with k_0 calculation");    

  

  // Spectrum options ------------------------------------------------------------------------------
  if (has_option("Emin"))
    _opt.Emin = get_option("Emin", 0.0);
  else
     throw InitFailedException("Optical Spectrum: Emin must be defined\n");

  if (has_option("Emax"))
    _opt.Emax = get_option("Emax", 0.0);
  else
    throw InitFailedException("Optical Spectrum: Emin must be defined\n");

  if (_opt.Emax < _opt.Emin)  throw InitFailedException("Optical Spectrum: Emax < Emin");

  if (has_option("dE"))
    _opt.dE = get_option("dE", 0.0);
  else
    throw InitFailedException("Optical Spectrum: dE must be defined\n");

  if (_opt.dE <= 0)  throw InitFailedException("Optical Spectrum: dE <= 0");


  _opt.get_occ = !get_option("compute_strengths",false);



  _opt.Gamma = get_option("broadening", 0.007);
  _opt.nr = get_option("refractive_index", 1.0);

  
  if (has_option("polarization") || has_option("polarisation"))
  {
    RealVectorValue polariz(3, 0.0);

    get_option("polarization", polariz);
    get_option("polarisation", polariz);

    _opt.polariz(1) = polariz(0);
    _opt.polariz(2) = polariz(1);
    _opt.polariz(3) = polariz(2);

    if (norm(_opt.polariz) != 0)
      _opt.polariz = _opt.polariz/norm(_opt.polariz);
    else
      throw InitFailedException("Optical Spectrum: polarization vector must be non zero");
  }

}


//=================================================================================

void Optics::do_init()
{
  //initial state----------------
  if (has_option("initial_state_model"))
  {
    std::string quantum_model;
    quantum_model = get_option("initial_state_model" , "");
    _initial_state_model = dynamic_cast<EigenvalueProblem*> (find_simulation(quantum_model));
    if (_initial_state_model == NULL)
      throw InitFailedException("Optics: invalid initial_state_model " + quantum_model);
  }
  else
    throw InitFailedException("Optics: initial_state_model must be defined\n");
  //--------------------------

  //final state----------------
  if (has_option("final_state_model"))
  {
    std::string quantum_model;
    quantum_model = get_option("final_state_model" , "");
    _final_state_model = dynamic_cast<EigenvalueProblem*> (find_simulation(quantum_model));
    if (_final_state_model == NULL)
      throw InitFailedException("Optics: invalid final_state_model " + quantum_model);

  }
  else
    throw InitFailedException("Optics: final_state_model must be defined\n");


  parse_options();

  _energy_mesh = new Mesh(1);

  unsigned int num_elem = (int)((_opt.Emax - _opt.Emin)/_opt.dE);

  MeshTools::Generation::build_cube (*_energy_mesh,
				     num_elem, 0, 0,
				     _opt.Emin, _opt.Emax,
				     0, 0,
				     0, 0,
				     EDGE2);


  _spectrum_x.resize(3*_energy_mesh->n_nodes(),0.0);
  _spectrum_y.resize(3*_energy_mesh->n_nodes(),0.0);
  _spectrum_z.resize(3*_energy_mesh->n_nodes(),0.0);


 // Kintegration options ----------------------------------------------------------

  init_k_space_integration();

}

//=================================================================================

void Optics::init_k_space_integration(void)
{

  // maybe this stuff should be taken from the intial/final state models?


   ModelOptions kopts;

   if (get_options().has_submodel("k_integration"))
   {
     ModelOptions::submodel_iterator it(get_options().submodels_begin("k_integration"));
     kopts = it->second;
   }
   else
     kopts.set_option("gamma_point_calculation", true);
   
   kopts.set_option("mesh_units", get_mesh_units());
   int k_dim = 3 - get_mesh().mesh_dimension();
   k_dim = kopts.get_option("k_space_dimension", k_dim);
   if (job == BULKMATREL)
     k_dim = 3;
   kopts.set_option("k_space_dimension", k_dim);

   kopts.set_option("verbose", SimulationOptions::verbose() );

   // these are the real space lattice vectors, in nm
   // why pi? Because then the default max k becomes 1 ( = 2*pi/(2*a) ), and
   // k_max can be interpreted in nm^-1
   RealVectorValue a(M_PI, 0, 0), b(0, M_PI, 0), c(0, 0, M_PI);

   // if there is an atomistic structure, we can take the lattice vectors from it
   // (they come in Angstrom!)
   if (_initial_state_model->get_atomistic_structure() != NULL)
   {
     _initial_state_model->get_atomistic_structure()->get_lattice_vectors(a, b, c);
     a *= 0.1;
     b *= 0.1;
     c *= 0.1;
   }

   switch (k_dim)
   {
     case 1:
       kopts.set_option("r1", c);
       break;

     case 2:
       kopts.set_option("r1", b);
       kopts.set_option("r2", c);
       break;

     case 3:
       kopts.set_option("r1", a);
       kopts.set_option("r2", b);
       kopts.set_option("r3", c);
       break;

     default:
       break;
   }


   Messages m;
   m.info("Setting up k-space integration");
   m.indent();

   _k_integration = KspaceIntegration::create(this, &Optics::calculate_for_k_point, kopts);

   if (_k_integration == NULL)
      throw InitFailedException("Could not create k-integration");

   _k_integration->init();
}

//==============================================//

void Optics::do_k_space_integration(void)
{
    // k-integration calls calculate_for_k_point()

    _k_integration->solve();
}

//==============================================//
void Optics::compute_matrix_elements()
{

  int verbose = SimulationOptions::verbose();

  if (verbose > 0)
    Messages::info("calculation of matrix elements for dipole optical transition...");

  set_states();

  if (job == BULKMATREL)  
    calculate_matrix_bulk(); 
  else
  {
    do_compute_matrix_elements();
  }



}




//==============================================//
void Optics::set_states()
{
  _i_states = _initial_state_model->get_solution();

  _initial_indices.clear();

  //-------------------------------------------------------------------------------------
  {
    std::vector<ID> inum = _initial_state_model->get_state_indices(_initial_state_particle);

    get_option("initial_eigenstates", _initial_indices);

    if (_initial_indices.size() > 0) 
    {
      if (_initial_indices.size() <= inum.size())
      {
        for (ID i = 0; i < _initial_indices.size(); i++)
        {
          //std::cout << "initial state indeces = " << temp[i] << std::endl;
          if(_initial_indices[i]<0 || _initial_indices[i]>inum.size()) continue;
          _initial_state_numbers[_initial_indices[i]] = inum[_initial_indices[i]];
        }
      }
    }
    else
    {
      //std::cout << "initial states defaults "<<std::endl;
      _initial_indices.resize(inum.size());
      for (ID i = 0; i < inum.size(); i++)
      {
        //std::cout << "initial state indeces = " << inum[i] << std::endl;
        _initial_indices[i] = i;
        _initial_state_numbers[i] = inum[i];
      }      
      
    }

    if (_initial_state_numbers.size() == 0)
      throw InitFailedException("Optics: no states available in "
                                + get_option("_initial_state_model", ""));
  }

  //-------------------------------------------------------------------------------------
  _f_states = _final_state_model->get_solution();  

  _final_indices.clear();

  {
    std::vector<unsigned int> fnum = _final_state_model->get_state_indices(_final_state_particle);

    //for (ID i = 0; i < fnum.size(); i++)
    // {
    //  std::cout << "final state indeces = " << fnum[i] << std::endl;
    // }

    get_option("final_eigenstates", _final_indices);

    if (_final_indices.size() > 0) 
    {
      if (_final_indices.size() <= fnum.size())
      {
        for (ID i = 0; i < _final_indices.size(); i++)
        {
          //std::cout << "final state indeces = " << temp[i] << std::endl;
          if(_final_indices[i]<0 || _final_indices[i]>fnum.size()) continue;
          _final_state_numbers[_final_indices[i]] = fnum[_final_indices[i]];
        }
      }
    }
    else
    {
      //std::cout << "initial states defaults "<<std::endl;
      _final_indices.resize(fnum.size());
      for (ID i = 0; i < fnum.size(); i++)
      {
        //std::cout << "final state indeces = " << fnum[i] << std::endl;
        _final_indices[i] = i;
        _final_state_numbers[i] = fnum[i];
        
      }      
      
    }


    if (_final_state_numbers.size() == 0)
      throw InitFailedException("Optics: no states available in " 
                               + get_option("_final_state_model", "") );
  }
  //-------------------------------------------------------------------------------------

}  
//==============================================//

void Optics::calculate_for_k_point(const Point& k_point,
                                   DofField& spectrum,
                                   double& integrated_quantity)
{

  short k_dim = _k_integration->get_k_space_dimension();
  if (plot_solution("EigenStates_k_0"))
    TiberCad::prepend_to_filename_suffix("k_0");
  else
  {
    ostringstream os;
    os << "k_(";
    switch (k_dim)
    {
      case 3:
        os << k_point(0) << ",";

      case 2:
        os << k_point(1) << ",";

      case 1:
        os << k_point(2);
        break;

      default:
        os << "0";
        break;
    }
    os << ")";
    TiberCad::prepend_to_filename_suffix(os.str());
  }

  if (!get_option("skip_solve", false))
  {
    _initial_state_model->solve_for_kpoint(k_point); //calculate eigenstates

    if( _initial_state_model != _final_state_model)
    {
      _final_state_model->solve_for_kpoint(k_point); //calculate eigenstates
    }
  }

  set_k_point(k_point);

  compute_matrix_elements(); //calculate matrix elements of P operator


  if ((plot_solution("EigenStates_k_0") && (k_point.size() < 1e-6)) ||
      plot_solution("EigenStates"))
  {
    _initial_state_model->plot();
    if( _initial_state_model != _final_state_model)
      _final_state_model->plot();

    plot_globaldata();
  }

  TiberCad::drop_first_filename_suffix();


  if (norm(_opt.polariz) == 0.0)
  {
    spectrum.resize(0);
    spectrum.reserve(9 * _energy_mesh->n_nodes());

    Tensor1 polariz; 
    polariz(1)=1.0;  polariz(2)=0.0; polariz(3)=0.0; 
    calculate_spectrum(*_energy_mesh, _opt.Gamma, polariz, _spectrum_z);
    spectrum.insert(spectrum.end(), _spectrum_z.begin(), _spectrum_z.end());

    polariz(1)=0.0;  polariz(2)=1.0; polariz(3)=0.0; 
    calculate_spectrum(*_energy_mesh, _opt.Gamma, polariz, _spectrum_z);
    spectrum.insert(spectrum.end(), _spectrum_z.begin(), _spectrum_z.end());

    polariz(1)=0.0;  polariz(2)=0.0; polariz(3)=1.0; 
    calculate_spectrum(*_energy_mesh, _opt.Gamma, polariz, _spectrum_z);
    spectrum.insert(spectrum.end(), _spectrum_z.begin(), _spectrum_z.end());
  }
  else
    calculate_spectrum(*_energy_mesh, _opt.Gamma, _opt.polariz,  spectrum);


  //for integrated quantity I take a sum of the map
  integrated_quantity = 0.0;

  for (unsigned int k=0; k < _energy_mesh->n_nodes(); k++)
    integrated_quantity += abs(spectrum[k]);

}

//================================================================//

void Optics::do_solve()
{
 
  // backup current solutions in initial and final state models
  ID init_sol_id = _initial_state_model->remember_current_solution();
  ID final_sol_id = INVALID_ID;
  if( _initial_state_model != _final_state_model)
    final_sol_id = _final_state_model->remember_current_solution();


  //if (norm(_opt.polariz) != 0)
  //{
  //  std::cout<<"(Optics) polarization: ( "<<_opt.polariz(1)<<" "
  //           <<_opt.polariz(2)<<" "<<_opt.polariz(3)<<")"<<std::endl;
  //}

  if (plot_solution("optical_spectrum")) 
  {
    do_k_space_integration();
    const DofField& spectra = _k_integration->get_solution();

    unsigned int n = 3*_energy_mesh->n_nodes();
    _spectrum_x.assign(spectra.begin(), spectra.begin() + n);

    if (norm(_opt.polariz) == 0.0)
    {
      _spectrum_y.assign(spectra.begin() + n, spectra.begin() + 2*n);
      _spectrum_z.assign(spectra.begin() + 2*n, spectra.begin() + 3*n);
    }
  }
  //---------------------------------------------------------------------------------
  
  if (plot_solution("optical_spectrum_k_0"))
  {

    double dummy;
    Point k_point;
    for(short i=0; i<3; i++) k_point(i) = _k_vector[i];

    calculate_for_k_point(k_point, _spectrum_x, dummy);

    if (norm(_opt.polariz) == 0.0)
    {
      unsigned int n = 3*_energy_mesh->n_nodes();
      _spectrum_y.assign(_spectrum_x.begin() + n, _spectrum_x.begin() + 2*n);
      _spectrum_z.assign(_spectrum_x.begin() + 2*n, _spectrum_x.begin() + 3*n);
      _spectrum_x.resize(n);
    }
  }

  _total_power = 0;
  _recombination = 0;

  // integrate the spectrum to obtain total emitted power and total recombination/generation rate

  /* I would maybe prefer to be able to loop over the elements:
  MeshBase::iterator it(_energy_mesh.elements_begin());
  MeshBase::iterator end(_energy_mesh.elements_end());
  for ( ; it != end; ++it)
  {
    const Elem* el = *it;
    double dE = el.size();
    _total_power +=
  }
  */

  // the last node
  size_t N = _energy_mesh->n_nodes() - 1;

  // use trapez formula
  for (size_t i = 1; i < N; i++)
  {
    double dP = _opt.dE * (_spectrum_x[i] + _spectrum_y[i] + _spectrum_z[i]);
    _total_power += dP;
    _recombination += dP / _energy_mesh->point(i)(0);
  }
  double dP0 = 0.5 * _opt.dE * (_spectrum_x[0] + _spectrum_y[0] + _spectrum_z[0]);
  double dPN = 0.5 * _opt.dE * (_spectrum_x[N] + _spectrum_y[N] + _spectrum_z[N]);
  _total_power += dP0 + dPN;
  _recombination += dP0 / _energy_mesh->point(0)(0) + dPN / _energy_mesh->point(N)(0);

  // k-space is assumed to be in 1/nm, and all output quantities
  // are given in 1/cm
  double area_dim_factor = 1.0;
  unsigned int kdim = _k_integration->get_k_space_dimension();
  int dim = (job == BULKMATREL) ? 0 : 3 - kdim;
  switch (dim)
  {
    case 0:
//      area_dim_factor *= Constants::bohr_radius * 1e2;
      area_dim_factor *= 1e-9 * 1e2;

    case 1:
//      area_dim_factor *= Constants::bohr_radius * 1e2;
      area_dim_factor *= 1e-9 * 1e2;

    case 2:
//      area_dim_factor *= Constants::bohr_radius * 1e2;
      area_dim_factor *= 1e-9 * 1e2;

    default:
      break;
  }

  // scale back to SI units
  _total_power *= Constants::elementary_charge / (Constants::atomic_time * area_dim_factor);
  _recombination /= (Constants::atomic_time * area_dim_factor);

  std::string type = get_option("type", "emission");
  ostringstream os;
  os << "Total " << ((type == "emission") ? "emitted" : "absorbed") << " power: "
      << _total_power << " " << get_solution_descriptor(OpticalPower).units() << "\n";
  os << "Total " << ((type == "emission") ? "recombination" : "generation") << " rate: "
      << _recombination << " " << get_solution_descriptor(Recombination).units();
  Messages::info(os.str());


  _initial_state_model->set_to_remembered_solution(init_sol_id);
  _initial_state_model->delete_remembered_solution(init_sol_id);
  if( _initial_state_model != _final_state_model)
  {
    _final_state_model->set_to_remembered_solution(final_sol_id);
    _final_state_model->delete_remembered_solution(final_sol_id);
  }
}

//=====================================================================================================


void Optics::do_calculate_spectrum(const Mesh& Energy, double Gamma,const Tensor1& polariz,
                                   DofField& spectrum )
{

  int verbose = SimulationOptions::verbose();

  if (verbose > 3)
  {
    ostringstream os;
    os << "calculation of optical spectrum, e = ("
        << polariz(1) << ", " <<  polariz(2) << ", " << polariz(3)
        << ")";
    Messages::info(os.str());
  }

  std::vector<double> fs_eigen_values;
  std::vector<double> is_eigen_values;

  std::vector<double> fs_occupations;
  std::vector<double> is_occupations;




  unsigned int n1 =  _initial_state_numbers.size();
  unsigned int n2 =  _final_state_numbers.size();


  _initial_state_model->get_eigenvalues(_initial_state_particle, is_eigen_values);

  _final_state_model->get_eigenvalues(_final_state_particle, fs_eigen_values);

  if (_opt.get_occ)
  {
    _initial_state_model->get_populations(_initial_state_particle, is_occupations);
    
    _final_state_model->get_populations(_final_state_particle, fs_occupations);
  }


  // loop on  eigenstates
  //cout << "Energy nodes: " << Energy.n_nodes()<<endl; 
  spectrum.clear();
  spectrum.reserve(3 * Energy.n_nodes());
  for (unsigned int el=0; el < 3 * Energy.n_nodes(); el++)
          spectrum.push_back(0.0);
  // try:
  // spectrum.resize(3 * Energy.n_nodes(), 0.0);


  for (unsigned i = 0; i < n1; i++)  // "upper" states
  {

    for (unsigned j = 0; j < n2; j++)  // "lower" states
    {

      double trans_energy =  abs(is_eigen_values[_initial_indices[i]]
	                       - fs_eigen_values[_final_indices[j]]);

      double f1 = 1.0;
      double f2 = 1.0;
      if (_opt.get_occ)
      {
        f1 = is_occupations[_initial_indices[i]];   // occupation for  electron
        
        f2 = fs_occupations[_final_indices[j]]; // occupation for  holes
      }


      Complex Me = _P_matrix[0][i][j] * polariz(1) +
	           _P_matrix[1][i][j] * polariz(2) +
	           _P_matrix[2][i][j] * polariz(3);


      double c = 1.0/Constants::fine_structure_constant;
      double omega = trans_energy/Constants::Hartree;

      //This is the right formula, as f1 is electron occupation probability and f2 is hole occupation probability.
      //Note that it differs from usual literature where usually f1 and f2 states initial state and final state
      //occupation probability, so it's related to electrons and it becomes f1*(1-f2)

      //spectrum[elem] += 1 / (2 * M_PI ) * (omega * omega) /(c*c*c)  * Lorenzian * abs (Me) * abs (Me) * f1 * f2;
      double strength = 2 * _opt.nr * (omega * omega) / (c*c*c) * abs (Me) * abs (Me);

      double power = strength * f1 * f2;
      double stimulated = strength * (f1 + f2 - 1);
      double gain = stimulated * M_PI * M_PI * c * c /
          (_opt.nr * _opt.nr * omega * omega * omega);

      //Note(alex): The original factor 1/(2*PI*PI) was changed to 1/(2*PI) (see below)
      //            and finally multiplied by 4 PI for angular integration (17/10/2011).
      //
      //According to Chuang's book the recombination rate should contain a pre-factor
      //(including 2 for spin sum and 2 for polarization and 4 Pi for angle integration)
      //
      //    nr^2 w^2           pi e^2      2       8 nr w e^2 c
      // ---------------- --------------- --- = ----------------------
      // pi^2 hbar c^2     nr c m^2 e0 w   V    4 pi e0 hbar (m^2 c^4)
      //
      // Extracting the prefactor m/hbar from the P-matrix, and multipling by a factor (hbar w) to get power emitted
      // we get the following prefactor (in which V has been removed to get total power emitted):
      //
      //     8 nr (hbar w)^2 e^2/(4 pi e0)
      //  = ------------------------------
      //             (hbar c)^3  hbar
      //
      //  Expressed in atomic units, hbar=1, e=1, 4 pi e0=1, c=1/fine_struct.
      //
      //  So the formula above, multiplied by 4*Pi for angle integration, and a factor of 4 for spin/pol deg
      //  agrees to Chuang's only if the factor is nr/(2*PI) rather than 1/(2*PI*PI).
      //
      //  note that  nr is  still missing !
      //
      //  2012-03-23, (Matthias): added nr taken from input file
      //
      //  2012-10-26, (Matthias): the second PI in the original factor (see Note Alex) was coming
      //                          from the Lorenzian. It is now included in the formula of the latter.



      unsigned int n_energy = Energy.n_nodes();
      for (unsigned int el = 0; el < n_energy; el++)
      {

        double En =  Energy.point(el)(0); //elem->centroid()(0);

        double Lorenzian =  0.5*Gamma/( ( trans_energy - En) *  ( trans_energy - En)
                            + (0.5*Gamma)*(0.5*Gamma)) / M_PI * Hartree;

        spectrum[el] += power * Lorenzian;
        spectrum[el + n_energy] += stimulated * Lorenzian;
        spectrum[el + 2*n_energy] += gain * Lorenzian;

     }
    }
  }
}

//============================================================================
void Optics::do_plot()
{
  // get  spectrum calculation  options from  opticsKP model  section
  // for  calculation of  spectrum for a single k-point

  if (plot_solution("matrix_elements"))
  {
    plot_globaldata();
  }


  int n_sets = (norm(_opt.polariz) == 0.0) ? 9 : 3;

  vector<double> results(_energy_mesh->n_nodes() * n_sets,0.0);
  
  vector<string> names(n_sets);
    
  string filename;
  string format = get_option("output_format", "grace");

  DataOutput data_output(*_energy_mesh, format);
  data_output.set_output_directory(get_output_directory());

  filename = get_name() + "_spectrum";

  if (plot_solution("optical_spectrum_k_0"))
  {
    filename += "_k_0";
  }
  
  filename += TiberCad::get_filename_suffix();
  


  string dimension;
  string gaindim;
  double area_dim_factor = 1;
  double gain_factor = 1;
  //double length_factor = Constants::bohr_radius * 1e2;
  double length_factor = 1e-9 * 1e2;

  short k_dim = (job == BULKMATREL) ? 3 : _k_integration->get_k_space_dimension();
  if (k_dim == 1)
  {
    dimension = "/cm";
    gaindim = "cm";
    area_dim_factor  = length_factor;
    gain_factor = length_factor;
  }
  else if (k_dim == 2)
  {
    dimension = "/cm^2";
    gaindim = "-";
    area_dim_factor  = length_factor * length_factor;
  }
  else if (k_dim == 3)
  {
    dimension = "/cm^3";
    gaindim = "1/cm";
    area_dim_factor  = length_factor *
        length_factor *
        length_factor;
    gain_factor = 1 / length_factor;
  }



  if (plot_solution("optical_spectrum") || plot_solution("optical_spectrum_k_0"))
  {

    if (n_sets == 9)
    {
      names[0] = "spontaneous_power_density_Px[W/eV" + dimension + "]";
      names[1] = "spontaneous_power_density_Py[W/eV" + dimension + "]";
      names[2] = "spontaneous_power_density_Pz[W/eV" + dimension + "]";
      names[3] = "stimulated_power_density_Px[W/eV" + dimension + "]";
      names[4] = "stimulated_power_density_Py[W/eV" + dimension + "]";
      names[5] = "stimulated_power_density_Pz[W/eV" + dimension + "]";
      names[6] = "gain_Px[" + gaindim + "]";
      names[7] = "gain_Py[" + gaindim + "]";
      names[8] = "gain_Pz[" + gaindim + "]";
    }
    else // n_sets == 3
    {
      names[0] = "spontaneous_power_density[W/eV" + dimension + "]";
      names[1] = "stimulated_power_density[W/eV" + dimension + "]";
      names[2] = "gain[" + gaindim + "]";
    }

    int point = 0;
    
    unsigned int n_el = _energy_mesh->n_nodes();
    for(unsigned int el = 0; el < n_el; el++)
    {
      double value;
      
      //--x - polarization
      value = _spectrum_x[el];
      value /= Constants::atomic_time; //[1/second/((bohr_radius)^kdim)]
      value *= Constants::elementary_charge; //[J/(eV*second)/((bohr_radius)^kdim)]
      value /= area_dim_factor ; //[J/(eV*s)/((cm)^kdim)]
      results[n_sets*point + 0] = value;

      if (n_sets == 9)
      {
        value = _spectrum_x[el + n_el];
        value /= Constants::atomic_time; //[1/second/((bohr_radius)^kdim)]
        value *= Constants::elementary_charge; //[J/(eV*second)/((bohr_radius)^kdim)]
        value /= area_dim_factor ; //[J/(eV*s)/((cm)^kdim)]
        results[n_sets*point + 3] = value;

        value = _spectrum_x[el + 2*n_el];
        value *= gain_factor;
        results[n_sets*point + 6] = value;
      }

      
      //--y - polarization
      value = _spectrum_y[el];
      value /= Constants::atomic_time; //[1/second/((bohr_radius)^kdim)]
      value *= Constants::elementary_charge; //[J/(eV*second)/((bohr_radius)^kdim)]
      value /= area_dim_factor ; //[J/(eV*s)/((cm)^kdim)]
      results[n_sets*point + 1] = value;
      
      if (n_sets == 9)
      {
        value = _spectrum_y[el + n_el];
        value /= Constants::atomic_time; //[1/second/((bohr_radius)^kdim)]
        value *= Constants::elementary_charge; //[J/(eV*second)/((bohr_radius)^kdim)]
        value /= area_dim_factor ; //[J/(eV*s)/((cm)^kdim)]
        results[n_sets*point + 4] = value;

        value = _spectrum_y[el + 2*n_el];
        value *= gain_factor;
        results[n_sets*point + 7] = value;
      }


      //--z - polarization
      value = _spectrum_z[el];
      value /= Constants::atomic_time; //[1/second/((bohr_radius)^kdim)]
      value *= Constants::elementary_charge; //[J/(eV*second)/((bohr_radius)^kdim)]
      value /= area_dim_factor ; //[J/(eV*s)/((cm)^kdim)]
      results[n_sets*point + 2] = value;

      if (n_sets == 9)
      {
        value = _spectrum_z[el + n_el];
        value /= Constants::atomic_time; //[1/second/((bohr_radius)^kdim)]
        value *= Constants::elementary_charge; //[J/(eV*second)/((bohr_radius)^kdim)]
        value /= area_dim_factor ; //[J/(eV*s)/((cm)^kdim)]
        results[n_sets*point + 5] = value;

        value = _spectrum_z[el + 2*n_el];
        value *= gain_factor;
        results[n_sets*point + 8] = value;
      }

      point++;
    }
       
    data_output.write_nodal_data(filename, results, names);

  }

}

//============================================================================

void Optics::print_info()
{

  int verbose = SimulationOptions::verbose();

  std::vector<double> is_pop, fs_pop;
  std::vector<double> is_ene, fs_ene;

  _initial_state_model->get_populations(_initial_state_particle, is_pop);
  _final_state_model->get_populations(_final_state_particle, fs_pop);
  _initial_state_model->get_eigenvalues(_initial_state_particle, is_ene);
  _final_state_model->get_eigenvalues(_final_state_particle, fs_ene);

  unsigned int n1 =  _initial_state_numbers.size();
  unsigned int n2 =  _final_state_numbers.size();

  // write down matrix elements for debugging
  if (verbose > 0)
  {

    std::cout<<"i-f      Ei       Ef    |Px|^2       fi    ff"
	     <<std::endl;
    for (int i = 0; i < n1; i++)
    {
      for (int j = 0; j < n2; j++)
      {

	unsigned int is = _initial_state_numbers[_initial_indices[i]];
	unsigned int fs = _final_state_numbers[_final_indices[j]];

	std::cout << is << "-" << fs << ": " << std::flush;
	std::cout << is_ene[is] << ":" << fs_ene[fs] << " " << std::flush;

	std::cout << std::norm(_P_matrix[0][i][j]) << "   ";
	std::cout << is_pop[is] << ", " << fs_pop[fs] << std::endl;

      }
    }

    std::cout<<"i-f      Ei       Ef    |Py|^2       fi    ff"
	     <<std::endl;
    for (int i = 0; i < n1; i++)
    {
      for (int j = 0; j < n2; j++)
      {

	unsigned int is = _initial_state_numbers[_initial_indices[i]];
	unsigned int fs = _final_state_numbers[_final_indices[j]];

	std::cout << is << "-" << fs << ": " << std::flush;
	std::cout << is_ene[is] << ":" << fs_ene[fs] << " " << std::flush;

	std::cout << std::norm(_P_matrix[1][i][j]) << "   ";
	std::cout << is_pop[is] << ", " << fs_pop[fs] << std::endl;

      }
    }

    std::cout<<"i-f      Ei       Ef    |Pz|^2       fi    ff"
	     <<std::endl;
    for (int i = 0; i < n1; i++)
    {
      for (int j = 0; j < n2; j++)
      {

	unsigned int is = _initial_state_numbers[_initial_indices[i]];
	unsigned int fs = _final_state_numbers[_final_indices[j]];

	std::cout << is << "-" << fs << ": " << std::flush;
	std::cout << is_ene[is] << ":" << fs_ene[fs] << " " << std::flush;

	std::cout << std::norm(_P_matrix[2][i][j]) << "   ";
	std::cout << is_pop[is] << ", " << fs_pop[fs] << std::endl;

      }
    }


  }


}

//========================================================================================
/*
void Optics::check_states()
{

  // Find maximum state index of initial states
  unsigned int i;
  unsigned int n_i =_initial_state_numbers.size();
  std::cerr<<"num initial states: "<< n_i<< std::endl;
  unsigned int n_i_st =_initial_state_model->get_num_states(_initial_state_particle);
  std::cerr<<"num computed:  "<< n_i_st<< std::endl;

  bool resized = false;

  for(i=n_i; i>0; --i)
  {
    if(_initial_state_numbers[i-1] > n_i_st-1)
    {
      _initial_indices.pop_back();
      resized = true;

    }
  }

  if(resized)
  {
    for(i=0; i<_initial_state_numbers.size(); i++)
    {
      std::cerr<<_initial_state_numbers[i]<<" ";
    }
    std::cerr<<std::endl;
  }

  // Find maximum state index of final states
  unsigned int n_f = _final_eigen_state_numbers.size();
  std::cerr<<"num final states: "<< n_f<< std::endl;
  unsigned int n_f_st =_final_state_model->get_num_states(_final_state_particle);
  std::cerr<<"num computed:  "<< n_f_st<< std::endl;

  resized = false;

  for(i=n_f; i>0; --i)
  {
    if(_final_eigen_state_numbers[i-1] > n_f_st-1)
    {
      _final_eigen_state_numbers.pop_back();
      resized = true;
    }
  }

  if(resized)
  {
    std::cerr<< "final states redefined: ";
    for(i=0; i<_final_eigen_state_numbers.size(); i++)
    {
      std::cerr<<_final_eigen_state_numbers[i]<<" ";
    }
    std::cerr<<std::endl;
  }



}
*/

void
Optics::plot_globaldata(void)
{
  string outdir = get_output_directory();

  string filename(outdir + "/" + get_output_filename() + ".dat");
  
  Messages::info("(OPT) write matrix elements in "+filename);

  ofstream file;
  file.open(filename.c_str());

  if (file.good())
  {
    file << "# initial_state final_state Px Py Pz\n";
    unsigned int n1 =  _initial_state_numbers.size();
    unsigned int n2 =  _final_state_numbers.size();

    for (unsigned int i = 0; i < n1; i++)
    {
      for (unsigned int j = 0; j < n2; j++)
      {
        file << i << "  " << j;
        for (unsigned int p = 0; p < 3; p++)
        {
          file << "  " << std::norm(_P_matrix[p][i][j]);
        }        
        file << "\n";
      }
    }
  }

}

