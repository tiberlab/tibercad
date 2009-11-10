// the following _HAS_ to be included first
#include "Read_MSH.h"

#include "ElectricalContact.h"
#include "Dopant.h"
#include "DriftDiffusion.h"
#include "DriftDiffusionProperties.h"

#include "Macrostrain.h"
#include "MacrostrainModelInterface.h"
#include "MacrostrainBoundaryProperties.h"

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

#include <algorithm>

using namespace std;
using namespace DriftDiffusionDefs;


void sweep_drain(double stop, int steps, DriftDiffusion& dd, double vg);


int main (int argc, char** argv)
{
  libMesh::init(argc, argv);
  {

    // general options
    GetPot input_file("options.in");

    const string searchpath = input_file("searchpath", ".");
    string meshfile = input_file("meshfile", "");
    string mesh_units = input_file("mesh_units", "1e-4");


    // drift-diffusion options
    GetPot dd_opt("dd.in");

    string curr_method = dd_opt("current_integration", "default");
    string dd_nonlin_rtol = dd_opt("nonlinear_tolerance", "1e-9");
    string dd_lin_rtol = dd_opt("linear_tolerance", "1e-12");
    string dd_integration_order = dd_opt("integration_order", "5");
    string dd_nonlin_max_it = dd_opt("nonlinear_max_it", "15");
    string dd_lin_max_it = dd_opt("linear_max_it", "500");
    string dd_ls_type = dd_opt("ls_type", "cubic");
    string dd_nonlin_ls_maxstep =
      dd_opt("nonlinear_ls_maxstep", "0.025");

    int characteristic = dd_opt("characteristic", 0);

    double vds_stop = dd_opt("vds_stop", 0.0);
    unsigned int vds_steps = dd_opt("vds_steps", 1);
    double vg_start = dd_opt("vg_start", 0.0);
    double vg_stop = dd_opt("vg_stop", 0.0);
    unsigned int vg_steps = dd_opt("vg_steps", 1);

    double vg = dd_opt("vg", 0.0);
    double vds = dd_opt("vds", 0.0);

    
    // strain options
    GetPot strain_opt("strain.in");
    unsigned int s_max_r_steps = strain_opt("max_r_steps", 4);
    int uniform_refinement = strain_opt("uniform_refinement", 0);
    double refine_fraction = strain_opt("refine_fraction", 0.8);
    double coarsen_fraction = strain_opt("coarsen_fraction", 0.0);
    unsigned int max_ref_level = strain_opt("max_ref_level", 10);
    double s_tolerance = strain_opt("tolerance", 1e-12);
    unsigned int max_shape_steps = strain_opt("max_shape_steps", 0);
    unsigned int substr_mat = strain_opt("substrate_material", 0);
    bool grown_on_substrate = strain_opt("grown_on_substrate", 0);
    double min_coord_substrate[3];
    min_coord_substrate[0] = strain_opt("xmin_s", 0.0);
    min_coord_substrate[1] = strain_opt("ymin_s", 0.0);
    min_coord_substrate[2] = strain_opt("zmin_s", 0.0);
    double max_coord_substrate[3];
    max_coord_substrate[0] = strain_opt("xmax_s", 0.0);
    max_coord_substrate[1] = strain_opt("ymax_s", 0.0);
    max_coord_substrate[2] = strain_opt("zmax_s", 0.0);
    unsigned int read_regions_from_mesh =
      input_file("read_regions_from_mesh", 1);
    bool calculate_atom_displacements =
      input_file("calculate_atom_displacements",false);
    std::string atom_structure_filename =
      input_file("atom_structure_filename","");
    std::string atom_displacements_filename =
      input_file("atom_displacements_filename","");

    
    unsigned int number_of_regions = 4;
    vector<unsigned int> phys_reg_ID(number_of_regions);
    phys_reg_ID[0] = 1; // AlGaN
    phys_reg_ID[1] = 2; // GaN (channel)
    phys_reg_ID[2] = 3; // AlGaN doped
    phys_reg_ID[3] = 4; // GaN doped
    
    vector<unsigned int> BC_reg_ID(4);
    BC_reg_ID[0] = 1; // gate
    BC_reg_ID[1] = 2; // source
    BC_reg_ID[2] = 3; // drain
    BC_reg_ID[3] = 4; // substrate

    
    Mesh mesh(2);
    MeshData_elements meshdata(mesh);
    meshdata.enable_compatibility_mode();
    
    cerr << "Read meshfile: " << meshfile << "\n";
    Read_MSH readmesh(meshfile, phys_reg_ID, BC_reg_ID, 2, mesh, meshdata);
    map<unsigned int, vector<unsigned int> > boundary_nodes;
    readmesh.get_BC_data(boundary_nodes);

    MeshUtils::assign_subdomain_ids(mesh, meshdata);
    
    mesh.print_info();

    // The device
    Device dev(mesh, boundary_nodes);

    // The materials
    Database d;
    d.set_search_path(searchpath);
    Material::set_database(d);

    // Strain simulation
    SimulationInterface* strainsim;
    {
      ModelOptions opts;
      opts["name"] = "strain";
      opts["substrate"] = "substrate";
      strainsim = SimulationInterface::create("macrostrain", opts);
    }
    
    // DD simulation
    SimulationInterface* ddsim;
    {
      ModelOptions dd_opts;
      dd_opts["name"] = "driftdiff";
      dd_opts["nonlin_max_it"] = "100";
      dd_opts["lin_max_it"] = dd_lin_max_it;
      dd_opts["nonlin_rel_tol"] = dd_nonlin_rtol;
      dd_opts["lin_rel_tol"] = dd_lin_rtol;
      dd_opts["integration_order"] = dd_integration_order;
      dd_opts["ls_type"] = dd_ls_type;
      dd_opts["ls_maxstep"] = dd_nonlin_ls_maxstep;
      dd_opts["mesh_units"] = mesh_units;
      dd_opts["pc_type"] = "composite";
      dd_opts["current_integration"] = curr_method;
      ddsim = SimulationInterface::create("drift-diffusion", dd_opts);
    }

    Material* mat_ptr;
    {
      MacrostrainModelInterface* strainmodel;
      DriftDiffusionProperties* ddmodel;

      ModelOptions matopts;
      matopts["structure"] = "wz";
      ModelOptions ddopts;
      ddopts["strain_simulation"] = "strain";
      ddopts["statistics"] = "FD";
      ModelOptions recomb_opts;
      recomb_opts["tau_n"] = "1e-9";
      recomb_opts["tau_p"] = "1e-9";
      ModelOptions mob_opts;


      // AlGan undoped, ID = 1
      matopts["x"] = "0.3";
      mat_ptr = Material::create("AlGaN", matopts);
      strainmodel = MacrostrainModelInterface::create("macrostrain");
      ddmodel = DriftDiffusionProperties::create("strained", ddopts);
      ddmodel->add_dopant(new Dopant(5e15, 0.025, 2, Dopant::N_TYPE));
      ddmodel->add_recombination_model("SRH", recomb_opts);
      mob_opts["mu0"] = "800";
      ddmodel->set_electron_mobility_model("constant", mob_opts);
      mob_opts["mu0"] = "200";
      ddmodel->set_hole_mobility_model("constant", mob_opts);
      mat_ptr->add_model(ddmodel, ddsim->get_id());
      mat_ptr->add_model(strainmodel, strainsim->get_id());
      dev.set_material(mat_ptr, 1);

      // Gan undoped, ID = 1
      mat_ptr = Material::create("GaN", matopts);
      strainmodel = MacrostrainModelInterface::create("macrostrain");
      ddmodel = DriftDiffusionProperties::create("strained", ddopts);
      ddmodel->add_dopant(new Dopant(5e15, 0.025, 2, Dopant::N_TYPE));
      ddmodel->add_recombination_model("SRH", recomb_opts);
      mob_opts["mu0"] = "800";
      ddmodel->set_electron_mobility_model("constant", mob_opts);
      mob_opts["mu0"] = "200";
      ddmodel->set_hole_mobility_model("constant", mob_opts);
      mat_ptr->add_model(ddmodel, ddsim->get_id());
      mat_ptr->add_model(strainmodel, strainsim->get_id());
      dev.set_material(mat_ptr, 2);

      // AlGan doped, ID = 1
      matopts["x"] = "0.3";
      mat_ptr = Material::create("AlGaN", matopts);
      strainmodel = MacrostrainModelInterface::create("macrostrain");
      ddmodel = DriftDiffusionProperties::create("strained", ddopts);
      ddmodel->add_dopant(new Dopant(5e19, 0.025, 2, Dopant::N_TYPE));
      recomb_opts["tau_n"] = "1e-12";
      recomb_opts["tau_p"] = "1e-12";
      ddmodel->add_recombination_model("SRH", recomb_opts);
      mob_opts["mu0"] = "200";
      ddmodel->set_electron_mobility_model("constant", mob_opts);
      mob_opts["mu0"] = "50";
      ddmodel->set_hole_mobility_model("constant", mob_opts);
      mat_ptr->add_model(ddmodel, ddsim->get_id());
      mat_ptr->add_model(strainmodel, strainsim->get_id());
      dev.set_material(mat_ptr, 3);

      // Gan doped, ID = 1
      mat_ptr = Material::create("GaN", matopts);
      strainmodel = MacrostrainModelInterface::create("macrostrain");
      ddmodel = DriftDiffusionProperties::create("strained", ddopts);
      ddmodel->add_dopant(new Dopant(5e19, 0.025, 2, Dopant::N_TYPE));
      recomb_opts["tau_n"] = "1e-12";
      recomb_opts["tau_p"] = "1e-12";
      ddmodel->add_recombination_model("SRH", recomb_opts);
      mob_opts["mu0"] = "800";
      ddmodel->set_electron_mobility_model("constant", mob_opts);
      mob_opts["mu0"] = "200";
      ddmodel->set_hole_mobility_model("constant", mob_opts);
      mat_ptr->add_model(ddmodel, ddsim->get_id());
      mat_ptr->add_model(strainmodel, strainsim->get_id());
      dev.set_material(mat_ptr, 4);

    }
    dev.init();

    set<ID> regions;
    regions.insert(1);
    regions.insert(2);
    regions.insert(3);
    regions.insert(4);
    SimulationEnvironment env(dev, regions);
    
    {
      ModelOptions ctopts;
      ElectricalContact* ct;
      Boundary* bd;

      bd = new Boundary("gate");
      ctopts["barrier_height"] = "1.2";
      ct = ElectricalContact::create("schottky", ctopts);
      bd->add_boundary_properties(ct, ddsim->get_id());
      env.add_boundary(bd, BC_reg_ID[0]);
      ctopts.clear();

      bd = new Boundary("source");
      ct= ElectricalContact::create("ohmic", ctopts);
      bd->add_boundary_properties(ct, ddsim->get_id());
      env.add_boundary(bd, BC_reg_ID[1]);

      bd = new Boundary("drain");
      ct= ElectricalContact::create("ohmic", ctopts);
      bd->add_boundary_properties(ct, ddsim->get_id());
      env.add_boundary(bd, BC_reg_ID[2]);

      bd = new Boundary("substrate");
      ctopts["material"] = "GaN";
      MacrostrainBoundaryProperties* ctsubstrate =
        MacrostrainBoundaryProperties::create("substrate", ctopts);
      bd->add_boundary_properties(ctsubstrate, strainsim->get_id());
      env.add_boundary(bd, BC_reg_ID[3]);
    }
    
    env.init();
    ddsim->set_environment(&env);
    strainsim->set_environment(&env);
    strainsim->init();

    env.prepare_for_solve();

    cout << "Solving strain... \n" << flush;
    strainsim->solve();
    Macrostrain& strain = *static_cast<Macrostrain*>(strainsim);
    strain.output_piezo("output/pizeo.gmv");
    strain.output_strain("output/strain.gmv");

    ddsim->init();

    DriftDiffusion& dd = *static_cast<DriftDiffusion*>(ddsim);
    cout << "Solving drift-diffusion... \n" << flush;

    dd.set_simulation_voltage("gate", 0.0);
    dd.set_simulation_voltage("source", 0.0);
    dd.set_simulation_voltage("drain", 0.0);

    dd.solve();
    const Scaling& sc = dd.get_scaling();
    cout << "Scaling parameters:\n";
    cout << "     phi0: " << sc.get_potential_scaling() << "\n";
    cout << "     x0  : " << sc.get_length_scaling() << "\n";
    cout << "     mu0 : " << sc.get_mobility_scaling() << "\n";
    cout << "     C0  : " << sc.get_density_scaling() << "\n";
    cout << "     t0  : " << sc.get_time_scaling() << "\n";
    cout << "     R0  : " << sc.get_recombination_scaling() << "\n\n";
    dd.remember_current_solution();
    vector<double> densities;
    vector<string> names;
    dd.build_densities(densities, names);
    GMVIO(dd.get_mesh()).write_nodal_data("output/eq_densities.gmv",
        densities, names);
    dd.build_band_edges(densities, names);
    GMVIO(dd.get_mesh()).write_nodal_data("output/eq_band_edges.gmv",
        densities, names);

    {
      ModelOptions dd_opts;
      dd_opts["nonlin_max_it"] = dd_nonlin_max_it;
      dd_opts["coupling"] = "full";
      ddsim->set_options(dd_opts);
    }

    if (characteristic)
    {
      cout << "calculate output characteristic...\n" << flush;
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
          ddsim->set_options(dd_opts);
          dd.set_simulation_voltage("gate", 0.0);
          dd.set_simulation_voltage("drain", 0.0);
          dd.set_electron_fermi_level(0.0);
          dd.set_hole_fermi_level(0.0);
          dd.set_electric_potential(0.0);
          dd.guess_equilibrium();
          dd.solve();
          dd.remember_current_solution();
          dd_opts["nonlin_max_it"] = dd_nonlin_max_it;
          dd_opts["coupling"] = "full";
          ddsim->set_options(dd_opts);
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
          << (*curr.find(dd.get_environment().get_boundary("gate"))).second
          << " "
          << (*curr.find(dd.get_environment().get_boundary("drain"))).second
          << " "
          << (*curr.find(dd.get_environment().get_boundary("source"))).second
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
            << (*curr.find(dd.get_environment().get_boundary("gate"))).second
            << " "
            << (*curr.find(dd.get_environment().get_boundary("drain"))).second
            << " "
            << (*curr.find(dd.get_environment().get_boundary("source"))).second
            << "\n" << flush;
          restart = false;
        }
        while (it != voltages.begin());
      }

    }

    delete ddsim;
    delete strainsim;


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
         << (*curr.find(dd.get_environment().get_boundary("gate"))).second
         << " "
         << (*curr.find(dd.get_environment().get_boundary("drain"))).second
         << " "
         << (*curr.find(dd.get_environment().get_boundary("source"))).second
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




