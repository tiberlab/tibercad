// the following _HAS_ to be included first
#include "Read_MSH.h"

#include "Dopant.h"
#include "RecombinationModelInterface.h"
#include "MobilityModelInterface.h"
#include "ElectricalContact.h"
#include "DriftDiffusion.h"
#include "StrainedSemiconductorModel.h"

#include "Macrostrain.h"
#include "WzDDsemiconductor.h"
#include "WzRotatedCrystal.h"
#include "WzPiezoelectricity.h"
#include "WzStiffness.h"

#include "Material.h"
#include "Device.h"
#include "Boundary.h"
#include "SimulationEnvironment.h"
#include "Database.h"
#include "MeshUtils.h"

#include "mesh_data_elements.h"

#include "mesh.h"
#include "mesh_modification.h"
#include "mesh_data.h"
#include "mesh_generation.h"
#include "elem.h"
#include "getpot.h"
#include "gnuplot_io.h"

#include <algorithm>

using namespace std;
using namespace DriftDiffusionDefs;



double alloy(double a, double b, double xa, double bowing = 0.0)
{
  return b + (a - b) * xa - bowing * xa *  (1 - xa);
};


int main (int argc, char** argv)
{

  libMesh::init(argc, argv);
  {

    // general options
    GetPot input_file("options.in");

    string searchpath = input_file("searchpath", ".");
    string meshfile = input_file("meshfile", "");
    string mesh_units = input_file("mesh_units", "1e-7");

    int fully = input_file("fully_coupled", 1);

    double in_content = input_file("in_content", 0.15);

    // drift-diffusion options
    GetPot dd_opt("dd.in");
    int characteristic = dd_opt("characteristic", 0);
    double n_doping = dd_opt("n_doping", 0.0);
    double p_doping = dd_opt("p_doping", 0.0);

    double start_voltage = dd_opt("vg_start", 0.0);
    double stop_voltage = dd_opt("vg_stop", 0.0);
    unsigned int voltage_steps = dd_opt("vg_steps", 1);

    string dd_nonlin_rtol = dd_opt("nonlinear_tolerance", "1e-9");
    string dd_nonlin_atol = dd_opt("nonlinear_abs_tolerance", "1e-12");
    string dd_lin_rtol = dd_opt("linear_tolerance", "1e-6");
    string dd_lin_atol = dd_opt("linear_abs_tolerance", "1e-9");
    string integration_order = dd_opt("integration_order", "5");
    string dd_nonlin_max_it = dd_opt("nonlinear_max_it", "15");
    string dd_lin_max_it = dd_opt("linear_max_it", "500");
    string ls_type = dd_opt("ls_type", "cubic");
    string dd_nonlin_ls_maxstep =
      dd_opt("nonlinear_ls_maxstep", "0.025");
    
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

    
    unsigned int number_of_regions = 3;
    vector<unsigned int> phys_reg_ID(number_of_regions);
    phys_reg_ID[0] = 1; // GaN (ev. n-doped)
    phys_reg_ID[1] = 2; // AlInN or InGaN (QW)
    phys_reg_ID[2] = 3; // GaN (ev. p-doped)
    
    vector<unsigned int> BC_reg_ID(3);
    BC_reg_ID[0] = 1; // anode
    BC_reg_ID[1] = 2; // cathode

    
    unsigned int dim = 1;
    Mesh mesh(1);
    MeshTools::Generation::build_line(mesh, 
        2060,
        -51.5, 51.5,
        EDGE2);


    MeshData_elements meshdata(mesh);
    meshdata.enable_compatibility_mode();
    
    map<unsigned int, vector<unsigned int> > boundary_nodes;
    boundary_nodes[1] = vector<unsigned int>(1, 2060);
    boundary_nodes[2] = vector<unsigned int>(1, 0);

    {
      map<const Elem*, vector<Number> > m;
      MeshBase::const_element_iterator it = mesh.elements_begin();
      const MeshBase::const_element_iterator end = mesh.elements_end();
      for ( ; it != end; ++it)
      {
        const Elem* elem = *it;

        double y = elem->centroid()(0);
        if (y < -1.5)
          m[elem] = vector<Number>(1, 1);
        else if (y < 1.5)
          m[elem] = vector<Number>(1, 2);
        else
          m[elem] = vector<Number>(1, 3);
      }
      meshdata.insert_elem_data(m);

    }


    MeshUtils::assign_subdomain_ids(mesh, meshdata);
 
    mesh.print_info();

    Device dev(mesh, boundary_nodes);

    /*****************************************************
     *
     * Setup of strain parameters
     *
     *****************************************************/

    bool periodicity[3];

    
    Macrostrain::strain_param a_strain;
    Macrostrain::strain_param b_strain;
    Macrostrain::strain_param c_strain;
    WzPiezoelectricity a_piezo;
    WzPiezoelectricity b_piezo;
    WzPiezoelectricity c_piezo;


    std::map<unsigned int, Piezoelectricity*> piezodata;
    piezodata[1] = &a_piezo;
    piezodata[2] = &b_piezo;
    piezodata[3] = &c_piezo;



    WzRotatedCrystal a_cryst;
    WzRotatedCrystal b_cryst;
    WzRotatedCrystal c_cryst;

    a_cryst.set_xyz_mil_direction("y", 1,  0, -1, 0) ;
    a_cryst.set_xyz_mil_direction("z", 1,  -2, 1, 0) ;
    a_cryst.set_xyz_mil_direction("x",  0,  0, 0,  1) ;
    b_cryst.set_xyz_mil_direction("y", 1,  0, -1, 0) ;
    b_cryst.set_xyz_mil_direction("z", 1,  -2, 1, 0) ;
    b_cryst.set_xyz_mil_direction("x",  0,  0, 0,  1) ;
    c_cryst.set_xyz_mil_direction("y", 1,  0, -1, 0) ;
    c_cryst.set_xyz_mil_direction("z", 1,  -2, 1, 0) ;
    c_cryst.set_xyz_mil_direction("x",  0,  0, 0,  1) ;

    a_cryst.set_lat_const(0.3189, 0.5185);
    b_cryst.set_lat_const(alloy(0.3545, 0.3189, in_content),
        alloy(0.5703, 0.5185, in_content));
    c_cryst.set_lat_const(0.3189, 0.5185);

    a_strain.crystal = &a_cryst;
    b_strain.crystal = &b_cryst;
    c_strain.crystal = &c_cryst;

    WzStiffness a_stiff;
    WzStiffness b_stiff;
    WzStiffness c_stiff;

    // GaN
    a_stiff.set_moduli(390.0, 145.0, 106.0, 398.0, 105.0);
    c_stiff.set_moduli(390.0, 145.0, 106.0, 398.0, 105.0);
    //
    a_piezo.set_moduli(1.27, -0.35, -0.3, -0.034); // C/m^2
    c_piezo.set_moduli(1.27, -0.35, -0.3, -0.034); // C/m^2

    // InGaN
    // first value: InN, second value: GaN
    b_stiff.set_moduli(
        alloy(223.0, 390.0, in_content),
        alloy(115.0, 145.0, in_content),
        alloy(92.0, 106.0, in_content),
        alloy(224.0, 398.0, in_content),
        alloy(48.0, 105.0, in_content));

    b_piezo.set_moduli(
        alloy(0.97, 1.27, in_content),
        alloy(-0.57, -0.35, in_content),
        alloy(-0.4, -0.3, in_content),
        alloy(-0.042, -0.034, in_content, -0.037));


    a_strain.C_tensor = &a_stiff;
    b_strain.C_tensor = &b_stiff;
    c_strain.C_tensor = &c_stiff;

    std::vector<int> x_dir(4);
    std::vector<int> y_dir(4);

    y_dir[0] =  1;  y_dir[1] = 0;  y_dir[2] = -1;  y_dir[3] =  0;
    x_dir[0] =  0;  x_dir[1] = 0;  x_dir[2] =  0;  x_dir[3] =  1;


    a_strain.crystal->calculate_lat_consts();
    b_strain.crystal->calculate_lat_consts();
    c_strain.crystal->calculate_lat_consts();
    a_strain.crystal->calculate_rot_matrix_miller(x_dir, y_dir); 
    b_strain.crystal->calculate_rot_matrix_miller(x_dir, y_dir);
    c_strain.crystal->calculate_rot_matrix_miller(x_dir, y_dir); 
    a_strain.C_tensor->rotate_to_calc_system(a_strain.crystal->RotMatrix);
    b_strain.C_tensor->rotate_to_calc_system(b_strain.crystal->RotMatrix);
    c_strain.C_tensor->rotate_to_calc_system(c_strain.crystal->RotMatrix);

    std::map<unsigned int, Macrostrain::strain_param> strain_params;
    strain_params.insert(pair<ID, Macrostrain::strain_param>(1, a_strain));
    strain_params.insert(pair<ID, Macrostrain::strain_param>(2, b_strain));
    strain_params.insert(pair<ID, Macrostrain::strain_param>(3, c_strain));


    Macrostrain::options opt;

    opt.intermediate_output = false;

    opt.max_r_steps        = s_max_r_steps;
    opt.uniform_refinement = uniform_refinement;
    opt.refine_fraction    = refine_fraction;
    opt.coarsen_fraction   = coarsen_fraction;
    opt.max_ref_level      = max_ref_level;
    opt.tolerance          = s_tolerance;
    opt.max_shape_steps    = max_shape_steps;
    opt.grown_on_substrate = grown_on_substrate;
    opt.substr_mat         = substr_mat;

    opt.calculate_atom_displacements = calculate_atom_displacements;
    opt.atom_structure_filename = atom_structure_filename;
    opt.atom_displacements_filename = atom_displacements_filename;

    opt.periodicity[0] = strain_opt("x-periodic", 0);
    opt.periodicity[1] = strain_opt("y-periodic", 0);
    opt.periodicity[2] = strain_opt("z-periodic", 0);

    vector<double> fp(3);

    fp[0] = fp[1] = fp[2] = 0.0;
    opt.fixed_point1 = fp;
    fp[0] = 100;
    opt.fixed_point2 = fp;
    fp[1] = 0.0; fp[0] = fp[2] = 1.0;
    opt.fixed_point3 = fp;

    double stress_value = strain_opt( "stress",0.0);
    map<unsigned int, double> stress_map;
    stress_map[1] = stress_value;


    Macrostrain strain_calculation(opt, dev.get_equation_systems(),
        "strainsys");

    strain_calculation.define_substrate_bc(1);
    strain_calculation.define_BC_map(boundary_nodes);
    strain_calculation.define_stress_value(stress_map);

    strain_calculation.assign_mesh_data(meshdata);
    strain_calculation.define_strain_parameters(strain_params);
    strain_calculation.define_piezo_moduli(piezodata);


    cout << "Solving strain... \n" << flush;
    strain_calculation.solve();
    strain_calculation.output_piezo("Piezo.gmv");

    /*****************************************************/


    /****************************************************************
     *
     * Drift diffusion definitions
     *
     ****************************************************************/

    SimulationInterface* dd_ptr;
    {
      ModelOptions dd_opts;
      //dd_opts["nonlin_max_it"] = "100";
      dd_opts["nonlin_max_it"] = dd_nonlin_max_it;
      dd_opts["lin_max_it"] = dd_lin_max_it;
      dd_opts["nonlin_rel_tol"] = dd_nonlin_rtol;
      dd_opts["nonlin_abs_tol"] = dd_nonlin_atol;
      dd_opts["lin_rel_tol"] = dd_lin_rtol;
      dd_opts["lin_abs_tol"] = dd_lin_atol;
      dd_opts["integration_order"] = integration_order;
      dd_opts["ls_type"] = ls_type;
      dd_opts["ls_maxstep"] = dd_nonlin_ls_maxstep;
      dd_opts["mesh_units"] = mesh_units;
      //dd_opts["pc_type"] = "composite";
      //dd_opts["pc_type"] = "jacobi";
      dd_ptr = SimulationInterface::create("drift-diffusion", dd_opts);
    }
    DriftDiffusion& dd = *dynamic_cast<DriftDiffusion*>(dd_ptr);


    DriftDiffusionProperties* nside = 
      DriftDiffusionProperties::create("strained");

    DriftDiffusionProperties* pside = 
      DriftDiffusionProperties::create("strained");

    DriftDiffusionProperties* well = 
      DriftDiffusionProperties::create("strained");

    {
      StrainedSemiconductorModel* sc =
        static_cast<StrainedSemiconductorModel*>(nside);
      sc->set_macrostrain(&strain_calculation);
      sc = static_cast<StrainedSemiconductorModel*>(pside);
      sc->set_macrostrain(&strain_calculation);
      sc = static_cast<StrainedSemiconductorModel*>(well);
      sc->set_macrostrain(&strain_calculation);
    }

    nside->set_statistics(TiberCad::FERMIDIRAC);
    pside->set_statistics(TiberCad::FERMIDIRAC);
    well->set_statistics(TiberCad::FERMIDIRAC);


    nside->add_dopant(new Dopant(n_doping, 0.025, 2, Dopant::N_TYPE));
    pside->add_dopant(new Dopant(p_doping, 0.17, 4, Dopant::P_TYPE));


    ModelOptions opts;
    opts["tau_n"] = "1e-10";
    opts["tau_p"] = "1e-10";
    nside->add_recombination_model("SRH", opts);
    pside->add_recombination_model("SRH", opts);
    opts["tau_n"] = "1e-7";
    opts["tau_p"] = "1e-7";
    well->add_recombination_model("SRH", opts);
    opts.clear();
    
    opts["C"] = "1e-10";
    nside->add_recombination_model("direct", opts);
    pside->add_recombination_model("direct", opts);
    well->add_recombination_model("direct", opts);
    opts.clear();

    opts["mu0"] = "200";
    nside->set_electron_mobility_model("constant", opts);
    pside->set_electron_mobility_model("constant", opts);
    well->set_electron_mobility_model("constant", opts);

    opts["mu0"] = "50";
    nside->set_hole_mobility_model("constant", opts);
    pside->set_hole_mobility_model("constant", opts);
    well->set_hole_mobility_model("constant", opts);
    opts.clear();

    Database d;
    d.set_search_path(searchpath);
    Material::set_database(d);

    opts["structure"] = "wz";
    Material* mat_n = Material::create("GaN", opts);
    mat_n->add_model(nside, dd.get_id());
    mat_n->init();

    Material* mat_p = Material::create("GaN", opts);
    mat_p->add_model(pside, dd.get_id());
    mat_p->init();

    opts["x"] = "0.14";
    Material* mat_w = Material::create("InGaN", opts);
    mat_w->add_model(well, dd.get_id());
    mat_w->init();

    dev.set_material(mat_n, 1);
    dev.set_material(mat_w, 2);
    dev.set_material(mat_p, 3);

    set<ID> regions;
    regions.insert(1);
    regions.insert(2);
    regions.insert(3);
    SimulationEnvironment dd_env(dev, regions);

    ModelOptions ctopts;
    ctopts["zero_field"] = "true";
    ctopts["zero_grad_fermi_e"] = "true";
    ctopts["zero_grad_fermi_h"] = "false";
    ElectricalContact* anode = ElectricalContact::create("ohmic", ctopts);
    ctopts["zero_grad_fermi_e"] = "false";
    ctopts["zero_grad_fermi_h"] = "true";
    ElectricalContact* cathode = ElectricalContact::create("ohmic", ctopts);

    Boundary* bd_anode = new Boundary("anode");
    Boundary* bd_cathode = new Boundary("cathode");
    bd_anode->add_boundary_properties(anode, dd.get_id());
    bd_cathode->add_boundary_properties(cathode, dd.get_id());
    dd_env.add_boundary(bd_anode, 1);
    dd_env.add_boundary(bd_cathode, 2);

    dd_env.init();
    dd.set_environment(&dd_env);

    


    mesh.print_info();
    
  

    dd.enable_mesh_refinement();

    /*****************************************************/

    

    dd.init();

    cout << "Solving drift-diffusion... \n" << flush;

    dd.set_simulation_voltage("cathode", 0.0);
    dd.set_simulation_voltage("anode", 0.0);
    
    cout << "Solving equilibrium... " << flush;
    dd.guess_equilibrium();
    try { dd.solve(); }
    catch (...) {}
    const Scaling& sc = dd.get_scaling();
    cout << "Scaling parameters:\n";
    cout << "     phi0: " << sc.get_potential_scaling() << "\n";
    cout << "     x0  : " << sc.get_length_scaling() << "\n";
    cout << "     mu0 : " << sc.get_mobility_scaling() << "\n";
    cout << "     C0  : " << sc.get_density_scaling() << "\n";
    cout << "     t0  : " << sc.get_time_scaling() << "\n";
    cout << "     R0  : " << sc.get_recombination_scaling() << "\n\n";

    dd.remember_current_solution();
    {
      vector<double> densities;
      vector<string> names;
      cout << "done (nr. iterations: " << dd.get_n_nonlinear_iterations() <<
        ", final residual: " << dd.get_final_residual() << ")\n" << flush;

      dd.build_densities(densities, names);
      GnuPlotIO(dd.get_mesh(), "Equilibrium densities",
          GnuPlotIO::GRID_ON).write_nodal_data("output/densities_eq",
            densities, names);
      dd.build_band_edges(densities, names);
      GnuPlotIO(dd.get_mesh(), "Equilibrium band diagram",
          GnuPlotIO::GRID_ON).write_nodal_data("output/bands_eq",
            densities, names);
    }


    cout << "\nBegin sweep...\n" << flush;
    {
      ModelOptions dd_opts;
      dd_opts["nonlin_max_it"] = dd_nonlin_max_it;
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
    int n = voltage_steps + 1;
    if (step < 1e-6) n = 1;
    vector<double> voltages(n);
    for (int i = 0; i <= voltage_steps; i++)
    {
      voltages[i] = start_voltage + i * step;
    }

    vector<double>::iterator first_positive =
      find_if(voltages.begin(), voltages.end(),
          bind2nd(greater<double>(), delta_v));

    ofstream file;
    file.open("output/iv_char.dat");
    file << "# V      A/cm^2\n" << flush;

    vector<double>::iterator it = first_positive;
    vector<double>::iterator zero =
      find_if(voltages.begin(), first_positive,
          bind2nd(greater<double>(), -delta_v));

    if (zero != first_positive)
      file << "0.0 0.0 0.0 0.0\n";


    for ( ; it != voltages.end(); ++it)
    {
      dd.set_simulation_voltage("anode", *it);
      cout << " Solving U = " << *it << " V ... " << flush;
      try { dd.solve(); }
      catch (...) {}
      cout << "done (nr. iterations: " << dd.get_n_nonlinear_iterations() <<
        ", final residual: " << dd.get_final_residual() << ")\n" << flush;

      vector<double> densities;
      vector<string> names;
      ostringstream filename_d;
      filename_d << "output/densities_" << *it;
      dd.build_densities(densities, names);
      GnuPlotIO(dd.get_mesh(), "Equilibrium densities",
          GnuPlotIO::GRID_ON).write_nodal_data(filename_d.str(),
            densities, names);

      ostringstream filename_b;
      filename_b << "output/bands_" << *it;
      dd.build_band_edges(densities, names);
      GnuPlotIO(dd.get_mesh(), "Equilibrium band diagram",
          GnuPlotIO::GRID_ON).write_nodal_data(filename_b.str(),
            densities, names);

      const map<const Boundary*, double>& curr =
        dd.get_boundary_currents();
      file << *it << "  "
           << (*curr.find(bd_cathode)).second << "  "
           << (*curr.find(bd_anode)).second << "\n" << flush;
      cerr << "    I = " << (*curr.find(bd_cathode)).second << " A/cm^2\n";
    }


    it = zero;
    if (it != voltages.begin())
    {
      dd.set_to_remembered_solution();
      do
      {
        --it;
        dd.set_simulation_voltage("anode", *it);
        cout << " Solving U = " << *it << " V ... " << flush;
        dd.solve();
        cout << "done (nr. iterations: " << dd.get_n_nonlinear_iterations() <<
          ", final residual: " << dd.get_final_residual() << ")\n" << flush;

        vector<double> densities;
        vector<string> names;
        ostringstream filename_d;
        filename_d << "output/densities_" << *it;
        dd.build_densities(densities, names);
        GnuPlotIO(dd.get_mesh(), "Equilibrium densities",
            GnuPlotIO::GRID_ON).write_nodal_data(filename_d.str(),
              densities, names);

        ostringstream filename_b;
        filename_b << "output/bands_" << *it;
        dd.build_band_edges(densities, names);
        GnuPlotIO(dd.get_mesh(), "Equilibrium band diagram",
            GnuPlotIO::GRID_ON).write_nodal_data(filename_b.str(),
              densities, names);

        const map<const Boundary*, double>& curr =
          dd.get_boundary_currents();
        file << *it << "  "
          << (*curr.find(bd_cathode)).second << "  "
          << (*curr.find(bd_anode)).second << "\n" << flush;
        cerr << "    I = " << (*curr.find(bd_cathode)).second << " A/cm^2\n";
      }
      while (it != voltages.begin());
    }
    file.close();
    
    delete dd_ptr;
  }

  return libMesh::close();
}



