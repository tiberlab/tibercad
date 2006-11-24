// the following _HAS_ to be included first
#include "Read_MSH.h"

#include "ElectricalContact.h"
#include "Dopant.h"
#include "DriftDiffusion.h"
#include "SemiconductorModel.h"
#include "RecombinationModelInterface.h"
#include "MobilityModelInterface.h"

#include "Material.h"
#include "Device.h"
#include "Boundary.h"
#include "SimulationEnvironment.h"
#include "Database.h"
#include "MeshUtils.h"

#include "mesh.h"
#include "mesh_modification.h"
#include "mesh_data.h"
#include "mesh_refinement.h"
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

    string searchpath = input_file("searchpath", ".");
    string meshfile = input_file("meshfile", "");
    string material = input_file("material", "Si");

    double temperature = input_file("temperature", 300.0);

    double start_voltage = input_file("start_voltage", 0.0);
    double stop_voltage = input_file("stop_voltage", 0.0);
    unsigned int voltage_steps = input_file("voltage_steps", 1);


    const string method = input_file("simulation_method", "NEWTON");
    const string statistics = input_file("statistics", "B");

    string nonlin_rtol = input_file("nonlinear_tolerance", "1e-9");
    string nonlin_atol = input_file("nonlinear_abs_tolerance", "1e-12");
    string nonlin_stol = input_file("nonlinear_step_tolerance", "1e-6");
    string lin_atol = input_file("linear_abs_tolerance", "1e-9");
    string lin_rtol = input_file("linear_tolerance", "1e-12");
    string integration_order = input_file("integration_order", "5");
    string nonlin_max_it = input_file("nonlinear_max_it", "15");
    string lin_max_it = input_file("linear_max_it", "500");
    string nonlin_ls_maxstep = input_file("nonlinear_ls_maxstep", "0.025");

    string mesh_units = input_file("mesh_units", "1e-4");

    double n_doping = input_file("n_doping", 1e18);
    double p_doping = input_file("p_doping", 1e18);


    unsigned int refinement_steps = input_file("max_refinement_steps", 0);
    double refine_frac = input_file("refine_fraction", 0.7);
    double coarsen_frac = input_file("coarsen_fraction", 0.3);

    vector<unsigned int> phys_reg_ID(5);
    phys_reg_ID[0] = 1; // n
    phys_reg_ID[1] = 2; // p
    phys_reg_ID[2] = 3; // i
    phys_reg_ID[3] = 4; // n+
    phys_reg_ID[4] = 5; // p+
    
    vector<unsigned int> BC_reg_ID(2);
    BC_reg_ID[0] = 1; // anode
    BC_reg_ID[1] = 2; // cathode

    unsigned int dim = 1;
    Mesh mesh(dim);
    MeshData_elements meshdata(mesh);
    meshdata.enable_compatibility_mode();
    
    Read_MSH readmesh(meshfile, phys_reg_ID, BC_reg_ID, dim, mesh, meshdata);
    map<unsigned int, vector<unsigned int> > boundary_nodes;
    readmesh.get_BC_data(boundary_nodes);
    

    MeshUtils::assign_subdomain_ids(mesh, meshdata);


    mesh.print_info();

    SimulationOptions::temperature = temperature;


    SimulationInterface* dd_ptr;
    {
      ModelOptions dd_opts;
      dd_opts["nonlin_max_it"] = "100";
      dd_opts["lin_max_it"] = lin_max_it;
      dd_opts["nonlin_rel_tol"] = nonlin_rtol;
      dd_opts["nonlin_step_tol"] = nonlin_stol;
      dd_opts["nonlin_abs_tol"] = nonlin_atol;
      dd_opts["lin_rel_tol"] = lin_rtol;
      dd_opts["integration_order"] = integration_order;
      dd_opts["ls_maxstep"] = nonlin_ls_maxstep;
      dd_opts["mesh_units"] = mesh_units;
      dd_opts["pc_type"] = "composite";
      dd_opts["ls_type"] = "quadratic";
      dd_ptr = SimulationInterface::create("drift-diffusion", dd_opts);
    }
    DriftDiffusion& dd = *dynamic_cast<DriftDiffusion*>(dd_ptr);


    DriftDiffusionProperties* nside = 
      DriftDiffusionProperties::create("unstrained");

    DriftDiffusionProperties* pside = 
      DriftDiffusionProperties::create("unstrained");

    DriftDiffusionProperties* intr = 
      DriftDiffusionProperties::create("unstrained");

    DriftDiffusionProperties* pplus = 
      DriftDiffusionProperties::create("unstrained");

    DriftDiffusionProperties* nplus = 
      DriftDiffusionProperties::create("unstrained");


    nside->add_dopant(new Dopant(n_doping, 0.025, 2, Dopant::N_TYPE));
    pside->add_dopant(new Dopant(p_doping, 0.01, 4, Dopant::P_TYPE));
    intr->add_dopant(new Dopant(1e15, 0.025, 2, Dopant::N_TYPE));
    nplus->add_dopant(new Dopant(1e20, 0.025, 2, Dopant::N_TYPE));
    pplus->add_dopant(new Dopant(1e20, 0.01, 4, Dopant::P_TYPE));

    ModelOptions opts;
    opts["tau_n"] = "1e-9";
    opts["tau_p"] = "3e-10";
    nside->add_recombination_model("SRH", opts);
    pside->add_recombination_model("SRH", opts);
    opts["tau_n"] = "1e-7";
    opts["tau_p"] = "3e-8";
    intr->add_recombination_model("SRH", opts);
    opts["tau_n"] = "1e-11";
    opts["tau_p"] = "3e-12";
    nplus->add_recombination_model("SRH", opts);
    pplus->add_recombination_model("SRH", opts);

    opts.clear();
    opts["mu0"] = "284";
    nside->set_electron_mobility_model("constant", opts);
    pside->set_electron_mobility_model("constant", opts);
    opts["mu0"] = "1400";
    intr->set_electron_mobility_model("constant", opts);
    opts["mu0"] = "130";
    nplus->set_electron_mobility_model("constant", opts);
    pplus->set_electron_mobility_model("constant", opts);

    opts["mu0"] = "160";
    nside->set_hole_mobility_model("constant", opts);
    pside->set_hole_mobility_model("constant", opts);
    opts["mu0"] = "400";
    intr->set_hole_mobility_model("constant", opts);
    opts["mu0"] = "90";
    nplus->set_hole_mobility_model("constant", opts);
    pplus->set_hole_mobility_model("constant", opts);


    Database d;
    d.set_search_path(searchpath);
    Material::set_database(d);

    Material* mat_n = Material::create(material);
    mat_n->add_model(nside, dd.get_id());
    mat_n->init();

    Material* mat_p = Material::create(material);
    mat_p->add_model(pside, dd.get_id());
    mat_p->init();

    Material* mat_i = Material::create(material);
    mat_i->add_model(intr, dd.get_id());
    mat_i->init();

    Material* mat_pplus = Material::create(material);
    mat_pplus->add_model(pplus, dd.get_id());
    mat_pplus->init();

    Material* mat_nplus = Material::create(material);
    mat_nplus->add_model(nplus, dd.get_id());
    mat_nplus->init();

    Device dev(mesh, boundary_nodes);
    dev.set_material(mat_n, 1);
    dev.set_material(mat_p, 2);
    dev.set_material(mat_i, 3);
    dev.set_material(mat_nplus, 4);
    dev.set_material(mat_pplus, 5);

    set<ID> regions;
    regions.insert(1);
    regions.insert(2);
    regions.insert(3);
    regions.insert(4);
    regions.insert(5);
    SimulationEnvironment dd_env(dev, regions);


    ModelOptions ctopts;
    ElectricalContact* anode = ElectricalContact::create("ohmic", ctopts);
    ElectricalContact* cathode = ElectricalContact::create("ohmic", ctopts);

    Boundary* bd_anode = new Boundary("anode");
    Boundary* bd_cathode = new Boundary("cathode");
    bd_anode->add_boundary_properties(anode, dd.get_id());
    bd_cathode->add_boundary_properties(cathode, dd.get_id());
    dd_env.add_boundary(bd_anode, 1);
    dd_env.add_boundary(bd_cathode, 2);


    dd_env.init();
    dd.set_environment(&dd_env);
    dd.init();


    dd.enable_mesh_refinement();

    
    dd.set_simulation_voltage("cathode", 0.0);
    dd.set_simulation_voltage("anode", 0.0);
    {
      vector<double> densities;
      vector<string> names;
      cout << "Solving equilibrium...\n" << flush;
      dd.guess_equilibrium();
      try { dd.solve(); }
      catch (...) {}
      dd.remember_current_solution();
      cout << "done (nr. iterations: " << dd.get_n_nonlinear_iterations() <<
        ", final residual: " << dd.get_final_residual() << ")\n" << flush;
      dd.build_densities(densities, names);
      GnuPlotIO(dd.get_mesh(), "Equilibrium densities",
          GnuPlotIO::GRID_ON).write_nodal_data("output/densities_eq",
          densities, names);
      dd.build_band_edges(densities, names);
      GnuPlotIO(dd.get_mesh(), "Equilibrium band profile",
          GnuPlotIO::GRID_ON).write_nodal_data("output/bands_eq",
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



    cout << "\nBegin sweep...\n" << flush;
    {
      ModelOptions dd_opts;
      dd_opts["nonlin_max_it"] = nonlin_max_it;
      dd_opts["coupling"] = "full";
      dd_ptr->set_options(dd_opts);
    }
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
      dd.enable_mesh_refinement();
      dd.solve();
      cout << "done (nr. iterations: " << dd.get_n_nonlinear_iterations() <<
        ", final residual: " << dd.get_final_residual() << ")\n" << flush;
      vector<double> densities;
      vector<string> names;
      ostringstream filename;
      filename << "output/bands_" << *it;
      ostringstream title;
      title << "Band profile for " << *it << " V";
      dd.build_band_edges(densities, names);
      GnuPlotIO(dd.get_mesh(), title.str(),
          GnuPlotIO::GRID_ON).write_nodal_data(filename.str(),
          densities, names);
      ostringstream filename_d;
      filename_d << "output/densities_" << *it;
      ostringstream title_d;
      title_d << "Densities for " << *it << " V";
      dd.build_densities(densities, names);
      GnuPlotIO(dd.get_mesh(), title_d.str(),
          GnuPlotIO::GRID_ON).write_nodal_data(filename_d.str(),
          densities, names);

      const map<const Boundary*, double>& curr =
        dd.get_boundary_currents();
      file << *it << "  "
           << (*curr.find(bd_cathode)).second << "  "
           << (*curr.find(bd_anode)).second << "\n" << flush;
      cerr << "    I = " << (*curr.find(bd_cathode)).second << " A/cm\n";
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
        ostringstream filename;
        filename << "output/bands_" << *it;
        ostringstream title;
        title << "Band profile for " << *it << " V";
        dd.build_band_edges(densities, names);
        GnuPlotIO(dd.get_mesh(), title.str(),
          GnuPlotIO::GRID_ON).write_nodal_data(filename.str(),
            densities, names);
        ostringstream filename_d;
        filename_d << "output/densities_" << *it;
        ostringstream title_d;
        title_d << "Densities for " << *it << " V";
        dd.build_densities(densities, names);
        GnuPlotIO(dd.get_mesh(), title_d.str(),
          GnuPlotIO::GRID_ON).write_nodal_data(filename_d.str(),
            densities, names);

        restart = false;
        const map<const Boundary*, double>& curr =
          dd.get_boundary_currents();
        file << *it << "  "
          << (*curr.find(bd_cathode)).second << "  "
          << (*curr.find(bd_anode)).second << "\n" << flush;
        cerr << "    I = " << (*curr.find(bd_cathode)).second << " A/cm\n";
      }
      while (it != voltages.begin());
    }
    
    file.close();

    delete dd_ptr;
  }

  return libMesh::close();
}


