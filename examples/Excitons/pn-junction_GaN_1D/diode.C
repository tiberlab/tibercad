// the following _HAS_ to be included first
#include "Read_MSH.h"

#include "OhmicContact.h"
#include "DriftDiffusion.h"
#include "Dopant.h"
#include "SemiconductorModel.h"
#include "ExcitonTransport.h"
#include "ExcitonModel.h"
#include "RecombinationModelInterface.h"
#include "ExcitonDissociation.h"

#include "Material.h"
#include "Device.h"
#include "Boundary.h"
#include "SimulationEnvironment.h"
#include "Database.h"
#include "MeshUtils.h"

#include "mesh.h"
#include "mesh_data.h"
#include "mesh_generation.h"
#include "elem.h"
#include "getpot.h"
#include "gnuplot_io.h"
#include "equation_systems.h"

#include <algorithm>

using namespace std;
using namespace DriftDiffusionDefs;




int main (int argc, char** argv)
{
  libMesh::init(argc, argv);
  {

    GetPot input_file("diode.in");

    string meshfile = input_file("meshfile", "");
    string material = input_file("material", "Si");
    string material_dir = input_file("material_dir", ".");
    material = material_dir + "/" + material + ".dat";

    const double mesh_units = input_file("mesh_units", 1e-4);
    int scaling = input_file("scaling", 1);
    bool linearize = input_file("linearize", false);

    double start_voltage = input_file("start_voltage", 0.0);
    double stop_voltage = input_file("stop_voltage", 0.0);
    unsigned int voltage_steps = input_file("voltage_steps", 1);

    double min_voltage_step = input_file("min_voltage_step", 1e-3);

    const string method = input_file("simulation_method", "NEWTON");
    const string statistics = input_file("statistics", "B");

    double nonlin_rtol = input_file("nonlinear_tolerance", 1e-9);
    double dd_nonlin_atol = input_file("nonlinear_abs_tolerance", 1e-12);
    double lin_rtol = input_file("linear_tolerance", 1e-12);
    double dd_lin_atol = input_file("linear_abs_tolerance", 1e-9);
    int integration_order = input_file("integration_order", 5);
    int nonlin_max_it = input_file("nonlinear_max_it", 15);
    int lin_max_it = input_file("linear_max_it", 500);
    double nonlin_ls_maxstep = input_file("nonlinear_ls_maxstep", 0.025);

    double n_doping = input_file("n_doping", 1e18);
    double p_doping = input_file("p_doping", 1e18);
    string mu_e = input_file("electron_mobility", "800");
    string mu_h = input_file("hole_mobility", "200");
    string tau_n = input_file("recombination_time_n", "1e-9");
    string tau_p = input_file("recombination_time_p", "1e-9");
    string C_direct = input_file("direct_recombination", "1e-8");
    string X_gen = input_file("exciton_generation", "1e-9");
    
    string trapping = input_file("trapping_probability", "0.5");
    double rec_rad = input_file("recombination_time_x_radiative", 1e-11);
    double rec_nonrad =
      input_file("recombination_time_x_nonradiative", 1e-11);
    
    double damping = input_file("damping", 1e-6);

    unsigned int refinement_steps = input_file("max_refinement_steps", 0);

    vector<unsigned int> phys_reg_ID(2);
    phys_reg_ID[0] = 1; // n
    phys_reg_ID[1] = 2; // p
    
    vector<unsigned int> BC_reg_ID(2);
    BC_reg_ID[0] = 1; // anode
    BC_reg_ID[1] = 2; // cathode

    
    int dim = 1;
    Mesh mesh(dim);
    MeshData_elements meshdata(mesh);
    meshdata.enable_compatibility_mode();
    
    Read_MSH readmesh(meshfile, phys_reg_ID, BC_reg_ID, dim, mesh, meshdata);
    map<unsigned int, vector<unsigned int> > boundary_nodes;
    readmesh.get_BC_data(boundary_nodes);

    MeshUtils::assign_subdomain_ids(mesh, meshdata);
    
    mesh.print_info();


    ModelOptions dd_opts;
    dd_opts["name"] = "driftdiff";
    SimulationInterface* dd_ptr =
      SimulationInterface::create("drift-diffusion", dd_opts);
    DriftDiffusion& dd = *dynamic_cast<DriftDiffusion*>(dd_ptr);

    ModelOptions ex_opts;
    ex_opts["name"] = "excitons";
    SimulationInterface* ex_ptr =
      SimulationInterface::create("exciton-transport", ex_opts);
    ExcitonTransport& ex = *dynamic_cast<ExcitonTransport*>(ex_ptr);

    DriftDiffusionProperties* nside = 
      DriftDiffusionProperties::create("unstrained");

    DriftDiffusionProperties* pside = 
      DriftDiffusionProperties::create("unstrained");


    if (statistics == "FD")
    {
      nside->set_statistics(TiberCad::FERMIDIRAC);
      pside->set_statistics(TiberCad::FERMIDIRAC);
    }

    nside->add_dopant(new Dopant(n_doping, 0.025, 2, Dopant::N_TYPE));
    pside->add_dopant(new Dopant(5e16, 0.025, 2, Dopant::N_TYPE));
    pside->add_dopant(new Dopant(p_doping, 0.17, 4, Dopant::P_TYPE));


    ModelOptions opts;
    opts["tau_n"] = tau_n;
    opts["tau_p"] = tau_n;
    nside->add_recombination_model("SRH", opts);
    opts["tau_n"] = tau_p;
    opts["tau_p"] = tau_p;
    pside->add_recombination_model("SRH", opts);
    
    opts.clear();
    opts["C"] = C_direct;
    nside->add_recombination_model("direct", opts);
    pside->add_recombination_model("direct", opts);
    
    opts.clear();
    opts["C"] = X_gen;
    nside->add_recombination_model("exciton_generation", opts);
    pside->add_recombination_model("exciton_generation", opts);

    opts.clear();
    opts["mu0"] = mu_e;
    nside->set_electron_mobility_model("constant", opts);
    pside->set_electron_mobility_model("constant", opts);

    opts["mu0"] = mu_h;
    nside->set_hole_mobility_model("constant", opts);
    pside->set_hole_mobility_model("constant", opts);


    Database d;
    Material::set_database(d);

    Material* mat_n = Material::create("Si");
    mat_n->add_model(nside, dd.get_id());
    mat_n->init();

    Material* mat_p = Material::create("Si");
    mat_p->add_model(pside, dd.get_id());
    mat_p->init();


    Device dev(mesh, boundary_nodes);
    dev.set_material(mat_n, 1);
    dev.set_material(mat_p, 2);

    set<ID> regions;
    regions.insert(1);
    regions.insert(2);
    SimulationEnvironment dd_env(dev, regions);
    SimulationEnvironment ex_env(dev, regions);



    OhmicContact anode;
    //SchottkyContact anode("anode");
    //anode.set_schottky_barrier_height(0.8);
    //anode.set_zero_derivative_bc(FERMIE);
    OhmicContact cathode;
    //cathode.set_zero_derivative_bc(FERMIH);

    Boundary bd_anode("anode"), bd_cathode("cathode");
    bd_anode.add_boundary_properties(&anode, dd.get_id());
    bd_cathode.add_boundary_properties(&cathode, dd.get_id());
    dd_env.add_boundary(&bd_anode, 1);
    dd_env.add_boundary(&bd_cathode, 2);

    dd_env.init();
    dd.set_environment(&dd_env);
    dd.init();
    
    //ex_env.init();
    //ex.set_environment(&dd_env);


    ExcitonModel& xpart = *ExcitonModel::create();
    xpart.set_exciton_generation_model("exciton_generation");
    // this is needed, but a momentary quirk
    xpart.set_recombination_times(rec_rad, rec_nonrad);
    xpart.set_binding_energy(20.4e-3);
    xpart.set_effective_mass(1.02);  // (*me) in exciton model
    xpart.set_mobility(1500);
    xpart.init();



  
    DriftDiffusion::Options& params = dd.get_options();
    params.max_refinement_steps = refinement_steps;
    params.solver_params.nonlinear_max_iterations = 100;
    params.solver_params.linear_max_iterations = lin_max_it;
    params.solver_params.nonlinear_tolerance = nonlin_rtol;
    params.solver_params.nonlinear_abs_tolerance = dd_nonlin_atol;
    params.solver_params.ls_maxstep = nonlin_ls_maxstep;
    params.solver_params.linear_tolerance = lin_rtol;
    params.solver_params.linear_abs_tolerance = dd_lin_atol;
    params.min_voltage_step = min_voltage_step;
    params.linearize_continuity_eq = linearize;
    params.integration_order =
      static_cast<libMeshEnums::Order>(integration_order);

    if (!scaling)
      params.scaling_type = Scaling::NONE;
    //else
    //  params.local_scaling = true;

    //params.artificial_drift = true;


    // mesh drawn in um
    params.mesh_units = mesh_units;

    //dd.enable_mesh_refinement();


    //ex.init();
    //ex.set_exciton_model(&xpart);
    //ex.set_initial_guess(1.5);

    ExcitonTransport::Options& exparams = ex.get_options();
    exparams.max_refinement_steps = refinement_steps;
    exparams.solver_params.nonlinear_max_iterations = 50;
    exparams.solver_params.linear_max_iterations = lin_max_it;
    exparams.solver_params.nonlinear_tolerance = 1e-19;
    exparams.solver_params.nonlinear_abs_tolerance = 1e-32;
    exparams.solver_params.ls_maxstep = nonlin_ls_maxstep;
    exparams.solver_params.linear_tolerance = lin_rtol;
    exparams.solver_params.linear_abs_tolerance = dd_lin_atol;
    exparams.solver_params.pc_type = PCCOMPOSITE;
    exparams.integration_order = 
      static_cast<libMeshEnums::Order>(5);

    exparams.mesh_units = mesh_units;

    params.solver_params.pc_type = PCILU;
    dd.set_simulation_voltage("cathode", 0.0);
    dd.set_simulation_voltage("anode", 0.0);
    {
      cout << "Solving equilibrium...\n" << flush;
      // is default:
      //params.coupling = POISSON;
      dd.guess_equilibrium();
      try { dd.solve(); }
      catch (...) {}

      //try { ex.solve(); }
      //catch (...) {}

      dd.remember_current_solution();
      cout << "done (nr. iterations: " << dd.get_n_nonlinear_iterations() <<
        ", final residual: " << dd.get_final_residual() << ")\n" << flush;

      vector<double> densities;
      vector<string> names;
      dd.build_densities(densities, names);
      GnuPlotIO(dd.get_mesh()).write_nodal_data("output/densities_eq",
            densities, names);
      dd.build_band_edges(densities, names);
      GnuPlotIO(dd.get_mesh()).write_nodal_data("output/bands_eq",
            densities, names);

      ex.build_densities(densities, names);
      GnuPlotIO(dd.get_mesh()).write_nodal_data("output/excitons_eq",
            densities, names);
    }

    const Scaling& sc = dd.get_scaling();
    cout << "Scaling parameters:\n";
    cout << "     phi0: " << sc.get_potential_scaling() << "\n";
    cout << "     x0  : " << sc.get_length_scaling() << "\n";
    cout << "     mu0 : " << sc.get_mobility_scaling() << "\n";
    cout << "     C0  : " << sc.get_density_scaling() << "\n";
    cout << "     t0  : " << sc.get_time_scaling() << "\n";
    cout << "     R0  : " << sc.get_recombination_scaling() << "\n\n";

    params.solver_params.nonlinear_max_iterations = nonlin_max_it;

    ostringstream s;
    s << damping;
    opts["damping"] = s.str();
    opts["exciton_simulation"] = "x_transport";
    //nside->add_recombination_model("exciton_dissociation", opts);
    //pside->add_recombination_model("exciton_dissociation", opts);

    cout << "\nBegin sweep...\n" << flush;
    params.solver_params.pc_type = PCCOMPOSITE;
    // make a voltage sweep
    double delta_v = 1e-6;

    if (stop_voltage < start_voltage)
    {
      double tmp = stop_voltage;
      stop_voltage = start_voltage;
      start_voltage = tmp;
    }
    // now stop_voltage > start_voltage

    double step = (stop_voltage - start_voltage) / voltage_steps;
    vector<double> voltages(voltage_steps + 1);
    for (int i = 0; i <= voltage_steps; i++)
    {
      voltages[i] = start_voltage + i * step;
    }

    vector<double>::iterator first_positive =
      find_if(voltages.begin(), voltages.end(),
          bind2nd(greater<double>(), delta_v));


    ofstream file;
    file.open("output/iv_char.dat");
    file << "# V      A/cm\n";
    

    vector<double>::iterator it = first_positive;
    for ( ; it != voltages.end(); ++it)
    {
      dd.set_simulation_voltage("anode", *it);
      cout << " Solving U = " << *it << " V ...\n" << flush;
      if (*it >= 1.5)
        dd.enable_mesh_refinement();
      
      params.coupling = FULLYCOUPLED;
      if (*it < 0.0)
      {
        dd.set_electron_fermi_level(0.0);
        dd.set_hole_fermi_level(*it);
        params.coupling = POISSON;
      }
      dd.solve();



      if (*it >= 0.5)
      {

        try { ex.solve(); }
        catch (...) {}

        dd.solve();
        try { ex.solve(); }
        catch (...) {}
        {
          damping = pow(damping, 0.25);
          ostringstream s;
          s << damping;
          opts["damping"] = s.str();
          //rma->set_model_options(opts);
          //rmb->set_model_options(opts);
        }
        dd.solve();
        try { ex.solve(); }
        catch (...) {}
        /*
        {
          damping = pow(damping, 0.25);
          ostringstream s;
          s << damping;
          opts["damping"] = s.str();
          rma->set_model_options(opts);
          rmb->set_model_options(opts);
        }
        dd.solve();
        try { ex.solve(); }
        catch (...) {}
        {
          damping = pow(damping, 0.25);
          ostringstream s;
          s << damping;
          opts["damping"] = s.str();
          rma->set_model_options(opts);
          rmb->set_model_options(opts);
        }
        dd.solve();
        try { ex.solve(); }
        catch (...) {}
        {
          damping = pow(damping, 0.25);
          ostringstream s;
          s << damping;
          opts["damping"] = s.str();
          rma->set_model_options(opts);
          rmb->set_model_options(opts);
        }
        */
      }
      

      cout << "done (nr. iterations: " << dd.get_n_nonlinear_iterations() <<
        ", final residual: " << dd.get_final_residual() << ")\n" << flush;
      ostringstream filename;
      vector<double> densities;
      vector<string> names;
      ostringstream f;
      f.precision(3);
      f << fixed << *it;
      dd.build_band_edges(densities, names);
      GnuPlotIO(dd.get_mesh()).write_nodal_data("output/bands_"+f.str(),
            densities, names);
      
      dd.build_densities(densities, names);
      GnuPlotIO(dd.get_mesh()).write_nodal_data("output/densities_"+f.str(),
            densities, names);

      ex.build_densities(densities, names);
      GnuPlotIO(dd.get_mesh()).write_nodal_data("output/excitons_"+f.str(),
            densities, names);
      

      const map<const Boundary*, double>& curr =
        dd.get_boundary_currents();
      file << *it << "  "
           << (*curr.find(&bd_cathode)).second << "  "
           << (*curr.find(&bd_anode)).second << "\n" << flush;
      cerr << "    I = " << (*curr.find(&bd_cathode)).second << " A/cm\n";
    }

    vector<double>::iterator zero =
      find_if(voltages.begin(), first_positive,
          bind2nd(greater<double>(), -delta_v));

    if (zero != first_positive)
      file << "0.0 0.0 0.0 0.0\n" << flush;

    it = zero;
    if (it != voltages.begin())
    {
      bool restart = true;
      do
      {
        --it;
        dd.set_simulation_voltage("anode", *it);
        cout << " Solving U = " << *it << " V ...\n" << flush;
        dd.enable_mesh_refinement();
        if (restart)
          dd.set_to_remembered_solution();
        dd.solve();
        cout << "done (nr. iterations: " << dd.get_n_nonlinear_iterations() <<
          ", final residual: " << dd.get_final_residual() << ")\n" << flush;
        vector<double> densities;
        vector<string> names;
        ostringstream f;
        f.precision(3);
        f << fixed << *it;
        dd.build_band_edges(densities, names);
        dd.build_band_edges(densities, names);
        GnuPlotIO(dd.get_mesh()).write_nodal_data("output/bands_"+f.str(),
            densities, names);
      
        dd.build_densities(densities, names);
        GnuPlotIO(dd.get_mesh()).write_nodal_data("output/densities_"+f.str(),
            densities, names);

        restart = false;

        const map<const Boundary*, double>& curr =
          dd.get_boundary_currents();
        file << *it << "  "
          << (*curr.find(&bd_cathode)).second << "  "
          << (*curr.find(&bd_anode)).second << "\n" << flush;
        cerr << "    I = " << (*curr.find(&bd_cathode)).second << " A/cm\n";
      }
      while (it != voltages.begin());
    }
    
    file.close();

  }

  return libMesh::close();
}


