// the following _HAS_ to be included first
#include "Read_MSH.h"
#include "ReadISEGrid.h"

#include "Dopant.h"
#include "DriftDiffusion.h"
#include "DriftDiffusionProperties.h"
#include "RecombinationModelInterface.h"
#include "ElectricalContact.h"

#include "Material.h"
#include "Device.h"
#include "Boundary.h"
#include "SimulationEnvironment.h"
#include "Database.h"
#include "MeshUtils.h"

#include "mesh.h"
#include "mesh_modification.h"
#include "mesh_data.h"
#include "mesh_generation.h"
#include "elem.h"
#include "getpot.h"
#include "gmv_io.h"
#include "tecplot_io.h"
#include "GMVIO_cell.h"

#include <algorithm>

using namespace std;
using namespace DriftDiffusionDefs;


void sweep_drain(double stop, int steps, DriftDiffusion& dd, double vg);

Boundary* bd_source;
Boundary* bd_drain;
Boundary* bd_gate;

int main (int argc, char** argv)
{
  libMesh::init(argc, argv);
  {

    GetPot input("mesfet.in");

    const string searchpath = input("searchpath", ".");
    string meshfile = input("meshfile", "");
    string format = input("format", "gmsh");
    string mesh_units = input("mesh_units", "1e-4");

    int fully = input("fully_coupled", 1);

    string schottky_barrier = input("schottky_barrier", "0.8");
    
    double vds_stop = input("vds_stop", 0.0);
    unsigned int vds_steps = input("vds_steps", 1);
    double vg_start = input("vg_start", 0.0);
    double vg_stop = input("vg_stop", 0.0);
    unsigned int vg_steps = input("vg_steps", 1);

    double vds = input("vds", 0.0);

    int characteristic = input("characteristic", 0);


    const string statistics = input("statistics", "B");


    string nonlin_rtol = input("nonlinear_tolerance", "1e-9");
    string lin_rtol = input("linear_tolerance", "1e-12");
    string integration_order = input("integration_order", "5");
    string nonlin_max_it = input("nonlinear_max_it", "15");
    string lin_max_it = input("linear_max_it", "500");
    string ls_type = input("ls_type", "cubic");
    string nonlin_ls_maxstep =
      input("nonlinear_ls_maxstep", "0.025");



    unsigned int refinement_steps =
      input("max_refinement_steps", 0);
    double refine_frac = input("refine_fraction", 0.7);
    double coarsen_frac = input("coarsen_fraction", 0.0);

    vector<unsigned int> phys_reg_ID(4);
    phys_reg_ID[0] = 1; // bulk
    phys_reg_ID[1] = 2; // channel
    phys_reg_ID[2] = 3; // contact1
    phys_reg_ID[3] = 4; // contact2
    
    vector<unsigned int> BC_reg_ID(3);
    BC_reg_ID[0] = 1; // source
    BC_reg_ID[1] = 2; // gate
    BC_reg_ID[2] = 3; // drain

    int dim = 2;
    Mesh mesh(dim);
    MeshData_elements meshdata(mesh);
    meshdata.enable_compatibility_mode();
    
    map<unsigned int, vector<unsigned int> > boundary_nodes;
    if (format == "ise")
    {
      ReadISEGrid readmesh(meshfile.c_str());
      mesh.read(readmesh.fname_xda,  &meshdata); 
      meshdata.read(readmesh.fname_xta);
      readmesh.get_BC_data(boundary_nodes);
    }
    else
    {
      Read_MSH readmesh(meshfile, phys_reg_ID, BC_reg_ID, dim, mesh, meshdata);
      readmesh.get_BC_data(boundary_nodes);
    }

    MeshUtils::assign_subdomain_ids(mesh, meshdata);

    mesh.print_info();



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
      //dd_opts["pc_type"] = "composite";
      dd_ptr = SimulationInterface::create("drift-diffusion", dd_opts);
    }


    
    DriftDiffusionProperties* sub = 
      DriftDiffusionProperties::create("unstrained");
    {
      sub->add_dopant(new Dopant(1e17, 0.01, 4, Dopant::P_TYPE));
      ModelOptions opts;
      opts["tau_n"] = "1e-6";
      opts["tau_p"] = "3e-7";
      sub->add_recombination_model("SRH", opts);
      opts.clear();
      opts["mu0"] = "560";
      sub->set_electron_mobility_model("constant", opts);
      opts["mu0"] = "280";
      sub->set_hole_mobility_model("constant", opts);
    }


    DriftDiffusionProperties* schottky = 
      DriftDiffusionProperties::create("unstrained");
    { 
      schottky->add_dopant(new Dopant(1e18, 0.025, 2, Dopant::N_TYPE));
      ModelOptions opts;
      opts["tau_n"] = "1e-7";
      opts["tau_p"] = "3e-8";
      schottky->add_recombination_model("SRH", opts);
      opts.clear();
      opts["mu0"] = "283.8";
      schottky->set_electron_mobility_model("constant", opts);
      opts["mu0"] = "160.3";
      schottky->set_hole_mobility_model("constant", opts);
    }


    DriftDiffusionProperties* contact = 
      DriftDiffusionProperties::create("unstrained");
    {
      contact->add_dopant(new Dopant(5e19, 0.025, 2, Dopant::N_TYPE));
      ModelOptions opts;
      opts["tau_n"] = "2e-9";
      opts["tau_p"] = "6e-10";
      contact->add_recombination_model("SRH", opts);
      opts.clear();
      opts["mu0"] = "70";
      contact->set_electron_mobility_model("constant", opts);
      opts["mu0"] = "54";
      contact->set_hole_mobility_model("constant", opts);
    }

    DriftDiffusionProperties* contact2 = 
      DriftDiffusionProperties::create("unstrained");
    {
      contact2->add_dopant(new Dopant(5e19, 0.025, 2, Dopant::N_TYPE));
      ModelOptions opts;
      opts["tau_n"] = "2e-9";
      opts["tau_p"] = "6e-10";
      contact2->add_recombination_model("SRH", opts);
      opts.clear();
      opts["mu0"] = "70";
      contact2->set_electron_mobility_model("constant", opts);
      opts["mu0"] = "54";
      contact2->set_hole_mobility_model("constant", opts);
    }


    Database d;
    d.set_search_path(searchpath);
    Material::set_database(d);


    ModelOptions opts;
    Material* mat_sub = Material::create("Si", opts);
    Material* mat_schottky = Material::create("Si", opts);
    Material* mat_contact = Material::create("Si", opts);
    Material* mat_contact2 = Material::create("Si", opts);

    mat_sub->add_model(sub, dd_ptr->get_id());
    mat_schottky->add_model(schottky, dd_ptr->get_id());
    mat_contact->add_model(contact, dd_ptr->get_id());
    mat_contact2->add_model(contact2, dd_ptr->get_id());

    mat_sub->init();
    mat_schottky->init();
    mat_contact->init();
    mat_contact2->init();

    Device dev(mesh, boundary_nodes);
    dev.set_material(mat_sub, 1);
    dev.set_material(mat_schottky, 2);
    dev.set_material(mat_contact, 3);
    dev.set_material(mat_contact2, 4);

    std::set<ID> ids;
    ids.insert(1);
    ids.insert(2);
    ids.insert(3);
    ids.insert(4);
    SimulationEnvironment dd_env(dev, ids);

    
    ModelOptions ctopts;
    ElectricalContact* source = ElectricalContact::create("ohmic", ctopts);
    ElectricalContact* drain = ElectricalContact::create("ohmic", ctopts);
    ctopts["schottky_barrier"] = schottky_barrier;
    ElectricalContact* gate = ElectricalContact::create("schottky", ctopts);

    bd_source = new Boundary("source");
    bd_drain = new Boundary("drain");
    bd_gate = new Boundary("gate");
    bd_source->add_boundary_properties(source, dd_ptr->get_id());
    bd_drain->add_boundary_properties(drain, dd_ptr->get_id());
    bd_gate->add_boundary_properties(gate, dd_ptr->get_id());
    dd_env.add_boundary(bd_source, 1);
    dd_env.add_boundary(bd_gate, 2);
    dd_env.add_boundary(bd_drain, 3);


    dd_env.init();
    dd_ptr->set_environment(&dd_env);
    dd_ptr->init();
  

    //dd.enable_mesh_refinement();
    
    DriftDiffusion& dd = static_cast<DriftDiffusion&>(*dd_ptr);

    dd.set_simulation_voltage("source", 0.0);
    dd.set_simulation_voltage("gate", 0.0);
    dd.set_simulation_voltage("drain", 0.0);


    dd.guess_equilibrium();
    dd.solve();
    dd.remember_current_solution();

    const Scaling& sc = dd.get_scaling();
    cout << "Scaling parameters:\n";
    cout << "     phi0: " << sc.get_potential_scaling() << "\n";
    cout << "     x0  : " << sc.get_length_scaling() << "\n";
    cout << "     mu0 : " << sc.get_mobility_scaling() << "\n";
    cout << "     C0  : " << sc.get_density_scaling() << "\n";
    cout << "     t0  : " << sc.get_time_scaling() << "\n";
    cout << "     R0  : " << sc.get_recombination_scaling() << "\n\n";


    vector<double> densities;
    vector<string> names;
    dd.build_band_edges(densities, names);
    TecplotIO(dd.get_mesh()).write_nodal_data("output/eq_bands.plt",
        densities, names);

    dd.build_densities(densities, names);
    TecplotIO(dd.get_mesh()).write_nodal_data("output/eq_densities.plt",
        densities, names);


    {
      ModelOptions dd_opts;
      dd_opts["nonlin_max_it"] = nonlin_max_it;
      dd_opts["coupling"] = "full";
      dd_ptr->set_options(dd_opts);
    }


    if (characteristic)
    {

      // make a voltage sweep
      double delta_v = 1e-6;

      if (vg_stop < vg_start)
      {
        double tmp = vg_stop;
        vg_stop = vg_start;
        vg_start = tmp;
      }
      // now vg_stop > vg_start

      double step = (vg_stop - vg_start) / vg_steps;
      int n = vg_steps + 1;
      if (step < 1e-6) n = 1;
      vector<double> voltages(n);
      for (int i = 0; i < n; i++)
      {
        voltages[i] = vg_start + i * step;
      }

      vector<double>::iterator first_positive =
        find_if(voltages.begin(), voltages.end(),
            bind2nd(greater<double>(), -delta_v));

      map<double, vector<double> > iv_char;

      vector<double>::iterator it = first_positive;
      for ( ; it != voltages.end(); ++it)
      {
        dd.set_simulation_voltage("gate", *it);
        sweep_drain(vds_stop, vds_steps, dd, *it);
        dd.set_to_remembered_solution();
      }

      vector<double>::iterator zero =
        find_if(voltages.begin(), first_positive,
            bind2nd(greater<double>(), -delta_v));

      //if (zero != first_positive)
      //  sweep_drain(vds_stop, vds_steps, dd, 0.0);

      it = zero;
      if (it != voltages.begin())
      {
        // recalculate equilibrium
        {
          ModelOptions dd_opts;
          dd_opts["nonlin_max_it"] = "100";
          dd_opts["coupling"] = "poisson";
          dd_ptr->set_options(dd_opts);
          dd.set_simulation_voltage("gate", 0.0);
          dd.set_simulation_voltage("drain", 0.0);
          dd.set_electron_fermi_level(0.0);
          dd.set_hole_fermi_level(0.0);
          dd.set_electric_potential(0.0);
          dd.guess_equilibrium();
          dd.solve();
          dd.remember_current_solution();
          dd_opts["nonlin_max_it"] = nonlin_max_it;
          dd_opts["coupling"] = "full";
          dd_ptr->set_options(dd_opts);
        }

        do
        {
          --it;
          dd.set_simulation_voltage("gate", *it);
          sweep_drain(vds_stop, vds_steps, dd, *it);
          dd.set_to_remembered_solution();
        }
        while (it != voltages.begin());
      }
    }
    else
    {
      // Id-vg characteristic
      dd.set_simulation_voltage("drain", vds);
      dd.solve();
      dd.remember_current_solution();

      ostringstream filename;
      filename.precision(3);
      filename << "output/id_vg_" << fixed << vds << "V.dat";
      ofstream file;
      file.open(filename.str().c_str());
      file << "# Id-Vg characteristic\n# Vg\tIg\tId\tIs\n";
      
      // make a voltage sweep
      double delta_v = 1e-6;

      if (vg_stop < vg_start)
      {
        double tmp = vg_stop;
        vg_stop = vg_start;
        vg_start = tmp;
      }
      // now vg_stop > vg_start

      double step = (vg_stop - vg_start) / vg_steps;
      vector<double> voltages(vg_steps + 1);
      for (int i = 0; i <= vg_steps; i++)
      {
        voltages[i] = vg_start + i * step;
      }

      vector<double>::iterator first_positive =
        find_if(voltages.begin(), voltages.end(),
            bind2nd(greater<double>(), -delta_v));

      map<double, vector<double> > iv_char;

      vector<double>::iterator it = first_positive;
      for ( ; it != voltages.end(); ++it)
      {
        dd.set_simulation_voltage("gate", *it);
        cout << "Vgs = " << *it << " Vds = " << vds << "\n" << flush;
        dd.solve();
        const map<const Boundary*, double>& curr =
          dd.get_boundary_currents();
        file << *it << " "
          << (*curr.find(bd_gate)).second
          << " "
          << (*curr.find(bd_drain)).second
          << " "
          << (*curr.find(bd_source)).second
          << "\n" << flush;
      }

      vector<double>::iterator zero =
        find_if(voltages.begin(), first_positive,
            bind2nd(greater<double>(), -delta_v));

      it = zero;
      if (it != voltages.begin())
      {
        bool restart = true;
        do
        {
          --it;
          dd.set_simulation_voltage("gate", *it);
          cout << "Vgs = " << *it << " Vds = " << vds << "\n" << flush;
          if (restart)
            dd.set_to_remembered_solution();
          dd.solve();
          const map<const Boundary*, double>& curr =
            dd.get_boundary_currents();
          file << *it << " "
            << (*curr.find(bd_gate)).second
            << " "
            << (*curr.find(bd_drain)).second
            << " "
            << (*curr.find(bd_source)).second
            << "\n" << flush;
          restart = false;
        }
        while (it != voltages.begin());
      }

    }
  }

  return libMesh::close();
}


void sweep_drain(double stop, int steps,
    DriftDiffusion& dd, double vg)
{

  ostringstream filename;
  filename.precision(3);
  filename << "output/ids_" << fixed << vg << "V.dat";
  ofstream file;
  file.open(filename.str().c_str());
  file << "# Id-Vds characteristic for Vg = " << vg << "V\n# Vds\tIg\tId\tIs\n";

  // make a voltage sweep
  double delta_v = 1e-6;

  double start = 0.0;

  if (stop < 0.0) stop = -stop;

  double step = (stop - start) / steps;
  vector<double> voltages(steps + 1);
  for (int i = 0; i <= steps; i++)
  {
    voltages[i] = start + i * step;
  }

  vector<double> data;
  vector<string> names;

  bool remember = true;
  vector<double>::iterator it = voltages.begin();
  for ( ; it != voltages.end(); ++it)
  {
    dd.set_simulation_voltage("drain", *it);
    cout << "Vgs = " << vg << " Vds = " << *it << "\n" << flush;
    dd.enable_mesh_refinement();
    dd.solve();
    if (remember)
    {
      dd.remember_current_solution();
    }
    remember = false;
    const map<const Boundary*, double>& curr =
      dd.get_boundary_currents();
    file << *it << " "
         << (*curr.find(bd_gate)).second
         << " "
         << (*curr.find(bd_drain)).second
         << " "
         << (*curr.find(bd_source)).second
         << "\n" << flush;

    ostringstream f;
    f.precision(3);
    f << "_" << fixed << vg << "V_" << fixed << *it << "V.gmv";
    dd.build_band_edges(data, names);
    GMVIO(dd.get_mesh()).write_nodal_data("output/bands"+f.str(),
        data, names);
    dd.build_densities(data, names);
    GMVIO(dd.get_mesh()).write_nodal_data("output/densities"+f.str(),
        data, names);
    dd.build_electric_field(data, names);
    GMVIO_cell(dd.get_mesh()).write_ascii_cell_data("output/field"+f.str(),
        data, names);
    dd.build_current_density(data, names);
    GMVIO_cell(dd.get_mesh()).write_ascii_cell_data("output/current"+f.str(),
        data, names);
  }
  
  file.close();
}




