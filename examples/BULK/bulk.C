// the following _HAS_ to be included first
#include "Read_MSH.h"

#include "ElectricalContact.h"
#include "FermiLevelPinning.h"
#include "SchottkyContact.h"
#include "Dopant.h"
#include "DriftDiffusion.h"
#include "DriftDiffusionProperties.h"

#include "Material.h"
#include "Device.h"
#include "Boundary.h"
#include "SimulationEnvironment.h"
#include "Database.h"
#include "MeshUtils.h"

#include "mesh.h"
#include "equation_systems.h"
#include "mesh_modification.h"
#include "mesh_data.h"
#include "mesh_generation.h"
#include "elem.h"
#include "getpot.h"
#include "GMVIO_cell.h"

#include <algorithm>

using namespace std;
using namespace DriftDiffusionDefs;



int main (int argc, char** argv)
{
  libMesh::init(argc, argv);
  {

    GetPot input("bulk.in");

    const string searchpath = input("searchpath", ".");
    const string meshfile = input("meshfile", "");
    const int dim = input("dimension", 3);

    double temperature = input("temperature", 300.0);

    double start_voltage = input("start_voltage", 0.0);
    double stop_voltage = input("stop_voltage", 0.0);
    unsigned int voltage_steps = input("voltage_steps", 1);

    const string material = input("material", "Si");
    const string structure = input("crystal_structure", "zb");
    const string molar_fraction = input("molar_fraction", "0.2");

    const string method = input("simulation_method", "NEWTON");
    const string statistics = input("statistics", "B");

    string nonlin_rtol = input("nonlinear_tolerance", "1e-9");
    string lin_rtol = input("linear_tolerance", "1e-12");
    string integration_order = input("integration_order", "5");
    string nonlin_max_it = input("nonlinear_max_it", "15");
    string lin_max_it = input("linear_max_it", "500");
    string ls_type = input("ls_type", "cubic");
    string nonlin_ls_maxstep =
      input("nonlinear_ls_maxstep", "0.025");

    double n_doping = input("n_doping", 1e15);
    double p_doping = input("p_doping", 0.0);
    double polarization_x = input("polarization_x", 0.0);
    double polarization_y = input("polarization_y", 0.0);
    
    string mesh_units = input("mesh_units", "1e-4");


    unsigned int refinement_steps =
      input("max_refinement_steps", 0);
    double refine_frac = input("refine_fraction", 0.7);
    double coarsen_frac = input("coarsen_fraction", 0.3);

    vector<unsigned int> phys_reg_ID(1);
    phys_reg_ID[0] = 1; 
    
    vector<unsigned int> BC_reg_ID(2);
    BC_reg_ID[0] = 1; // anode
    BC_reg_ID[1] = 2; // cathode

    
    Mesh mesh(dim);
    MeshData_elements meshdata(mesh);
    meshdata.enable_compatibility_mode();
    
    cerr << "Read meshfile: " << meshfile << "\n";
    Read_MSH readmesh(meshfile, phys_reg_ID, BC_reg_ID, dim, mesh, meshdata);
    map<unsigned int, vector<unsigned int> > boundary_nodes;
    readmesh.get_BC_data(boundary_nodes);

    MeshUtils::assign_subdomain_ids(mesh, meshdata);

    mesh.print_info();
    SimulationOptions::temperature = temperature;

    SimulationInterface* dd_ptr;
    {
      ModelOptions dd_opts;
      dd_opts["name"] = "driftdiff";
      dd_opts["nonlin_max_it"] = "100";
      dd_opts["lin_max_it"] = lin_max_it;
      dd_opts["nonlin_rel_tol"] = nonlin_rtol;
      dd_opts["lin_rel_tol"] = lin_rtol;
      dd_opts["integration_order"] = integration_order;
      dd_opts["ls_type"] = ls_type;
      dd_opts["ls_maxstep"] = nonlin_ls_maxstep;
      dd_opts["mesh_units"] = mesh_units;
      dd_opts["pc_type"] = "composite";
      dd_ptr = SimulationInterface::create("drift-diffusion", dd_opts);
    }



    DriftDiffusionProperties* nside = 
      DriftDiffusionProperties::create("unstrained");

    if (statistics == "FD")
      nside->set_statistics(TiberCad::FERMIDIRAC);
    else
      nside->set_statistics(TiberCad::BOLTZMANN);


    nside->add_dopant(new Dopant(n_doping, 0.025, 2, Dopant::N_TYPE));
    nside->add_dopant(new Dopant(p_doping, 0.01, 4, Dopant::P_TYPE));


    ModelOptions opts;
    opts["tau_n"] = "1e-7";
    opts["tau_p"] = "3e-8";
    nside->add_recombination_model("SRH", opts);
    opts["mu0"] = "800";
    nside->set_electron_mobility_model("constant", opts);
    opts["mu0"] = "200";
    nside->set_hole_mobility_model("constant", opts);


    Database d;
    d.set_search_path(searchpath);
    Material::set_database(d);

    opts.clear();
    opts["x"] = molar_fraction;
    opts["structure"] = structure;
    Material* mat_ptr = Material::create(material, opts);
    Material& mat = *mat_ptr;
    
    // test copy of models
    //DriftDiffusionProperties* ddprop =
    //  dynamic_cast<DriftDiffusionProperties*>(nside->copy());
    //delete nside;
    //mat.add_model(ddprop, dd.get_id());
      
    mat.add_model(nside, dd_ptr->get_id());
    mat.init();

    Device dev(mesh, boundary_nodes);
    dev.set_material(&mat, 1);

    SimulationEnvironment dd_env(dev, 1);


    ModelOptions ctopts;
    ElectricalContact* anode = ElectricalContact::create("ohmic", ctopts);
    ElectricalContact* cathode = ElectricalContact::create("ohmic", ctopts);

    Boundary* bd_anode = new Boundary("anode");
    Boundary* bd_cathode = new Boundary("cathode");
    bd_anode->add_boundary_properties(anode, dd_ptr->get_id());
    bd_cathode->add_boundary_properties(cathode, dd_ptr->get_id());
    dd_env.add_boundary(bd_anode, 1);
    dd_env.add_boundary(bd_cathode, 2);


    dd_env.init();
    dd_ptr->set_environment(&dd_env);
    dd_ptr->init();

    
    DriftDiffusion& dd = static_cast<DriftDiffusion&>(*dd_ptr);

    dd.set_simulation_voltage("cathode", 0.0);
    dd.set_simulation_voltage("anode", 0.0);
    {
      vector<double> densities;
      vector<string> names;
      cerr << "Solving equilibrium... " << flush;
      dd.guess_equilibrium();
      try { dd.solve(); }
      catch (...) {}
      dd.remember_current_solution();
      cerr << "done (nr. iterations: " << dd.get_n_nonlinear_iterations() <<
        ", final residual: " << dd.get_final_residual() << ")\n" << flush;
      dd.build_band_edges(densities, names);
      GMVIO(dd.get_mesh()).write_nodal_data("output/bands_eq.gmv",
          densities, names);
      dd.build_densities(densities, names);
      GMVIO(dd.get_mesh()).write_nodal_data("output/densities_eq.gmv",
          densities, names);
    }

    const Scaling& sc = dd.get_scaling();
    cerr << "Scaling parameters:\n";
    cout << "     phi0: " << sc.get_potential_scaling() << "\n";
    cout << "     x0  : " << sc.get_length_scaling() << "\n";
    cout << "     mu0 : " << sc.get_mobility_scaling() << "\n";
    cout << "     C0  : " << sc.get_density_scaling() << "\n";
    cout << "     t0  : " << sc.get_time_scaling() << "\n";
    cout << "     R0  : " << sc.get_recombination_scaling() << "\n\n";



    cout << "\nBegin sweep...\n" << flush;

    // calculate for electron and holes
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
      //if (*it >= 1.7)
      //  dd.enable_mesh_refinement();
      
      dd.solve();
      
      cout << "done (nr. iterations: " << dd.get_n_nonlinear_iterations() <<
        ", final residual: " << dd.get_final_residual() << ")\n" << flush;
      ostringstream filename;
      vector<double> densities;
      vector<string> names;
      ostringstream filename_b;
      filename_b << "output/bands_" << *it << ".gmv";
      dd.build_band_edges(densities, names);
      GMVIO(dd.get_mesh()).write_nodal_data(filename_b.str(),
          densities, names);
      ostringstream filename_d;
      filename_d << "output/densities_" << *it << ".gmv";
      dd.build_densities(densities, names);
      GMVIO(dd.get_mesh()).write_nodal_data(filename_d.str(),
          densities, names);
      ostringstream filename_f;
      filename_f << "output/field_" << *it << ".gmv";
      dd.build_electric_field(densities, names);
      GMVIO_cell(dd.get_mesh()).write_ascii_cell_data(filename_f.str(),
          densities, names);
      ostringstream filename_c;
      filename_c << "output/current_" << *it << ".gmv";
      dd.build_current_density(densities, names);
      GMVIO_cell(dd.get_mesh()).write_ascii_cell_data(filename_c.str(),
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
      dd.set_to_remembered_solution();
      do
      {
        --it;
        dd.set_simulation_voltage("anode", *it);
        cout << " Solving U = " << *it << " V ...\n" << flush;
        dd.enable_mesh_refinement();
        dd.solve();
        cout << "done (nr. iterations: " << dd.get_n_nonlinear_iterations() <<
          ", final residual: " << dd.get_final_residual() << ")\n" << flush;
        vector<double> densities;
        vector<string> names;
        ostringstream filename_b;
        filename_b << "output/bands_" << *it << ".gmv";
        dd.build_band_edges(densities, names);
        GMVIO(dd.get_mesh()).write_nodal_data(filename_b.str(),
            densities, names);
        ostringstream filename_d;
        filename_d << "output/densities_" << *it << ".gmv";
        dd.build_densities(densities, names);
        GMVIO(dd.get_mesh()).write_nodal_data(filename_d.str(),
            densities, names);
        ostringstream filename_f;
        filename_f << "output/field_" << *it << ".gmv";
        dd.build_electric_field(densities, names);
        GMVIO_cell(dd.get_mesh()).write_ascii_cell_data(filename_f.str(),
            densities, names);
        ostringstream filename_c;
        filename_c << "output/current_" << *it << ".gmv";
        dd.build_current_density(densities, names);
        GMVIO_cell(dd.get_mesh()).write_ascii_cell_data(filename_c.str(),
            densities, names);

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

