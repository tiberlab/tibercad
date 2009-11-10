// the following _HAS_ to be included first
#include "Read_MSH.h"

#include "ElementData.h"
#include "OhmicContact.h"
#include "BoundaryData.h"
#include "DDevice.h"
#include "DriftDiffusion.h"
#include "StrainedSemiconductorModel.h"

#include "macrostrain.h"
#include "mesh_data_elements.h"

#include "mesh.h"
#include "mesh_modification.h"
#include "mesh_data.h"
#include "mesh_generation.h"
#include "elem.h"
#include "getpot.h"
#include "gnuplot_io.h"
#include "WzDDsemiconductor.h"

#include <algorithm>

using namespace std;
using namespace DriftDiffusionDefs;


void set_boundary(BoundaryData& data, const vector<unsigned int>& nodes,
    ElectricalContact* desc, const Mesh& mesh);


double alloy(double a, double b, double xa, double bowing = 0.0)
{
  return b + (a - b) * xa - bowing * xa *  (1 - xa);
};

class Dummy {};

int main (int argc, char** argv)
{

  libMesh::init(argc, argv);
  {

    // general options
    GetPot input_file("options.in");

    string meshfile = input_file("meshfile", "");
    const double mesh_units = input_file("mesh_units", 1e-7);

    int fully = input_file("fully_coupled", 1);

    const string approx_order =
      input_file("approximation_order", "FIRST");

    double al_content = input_file("al_content", 0.48);
    double ga_content = input_file("ga_content", 0.47);
    int growth_dir = input_file("growth_dir",1);
    int solve_strain = input_file("solve_strain", 0);

    // drift-diffusion options
    GetPot dd_opt("dd.in");
    int characteristic = dd_opt("characteristic", 0);
    double n_doping = dd_opt("n_doping", 0.0);
    double p_doping = dd_opt("p_doping", 0.0);

    double start_voltage = dd_opt("vg_start", 0.0);
    double stop_voltage = dd_opt("vg_stop", 0.0);
    unsigned int voltage_steps = dd_opt("vg_steps", 1);

    const string method = dd_opt("simulation_method", "NEWTON");
    const string statistics = dd_opt("statistics", "FD");
    double dd_nonlin_rtol = dd_opt("nonlinear_tolerance", 1e-9);
    double dd_nonlin_atol = dd_opt("nonlinear_abs_tolerance", 1e-12);
    double dd_lin_rtol = dd_opt("linear_tolerance", 1e-6);
    double dd_lin_atol = dd_opt("linear_abs_tolerance", 1e-9);
    int integration_order = dd_opt("integration_order", 5);
    int dd_nonlin_max_it = dd_opt("nonlinear_max_it", 15);
    int dd_lin_max_it = dd_opt("linear_max_it", 500);
    double dd_nonlin_ls_maxstep =
      dd_opt("nonlinear_ls_maxstep", 0.025);
    unsigned int dd_max_r_steps =
      dd_opt("max_refinement_steps", 0);
    
    // strain options
    GetPot strain_opt("strain.in");
    unsigned int s_max_r_steps = strain_opt("max_r_steps", 4);
    int uniform_refinement = strain_opt("uniform_refinement", 0);
    double refine_fraction = strain_opt("refine_fraction", 0.8);
    double coarsen_fraction = strain_opt("coarsen_fraction", 0.0);
    unsigned int max_ref_level = strain_opt("max_ref_level", 10);
    double s_tolerance = strain_opt("tolerance", 1e-12);
    unsigned int max_shape_steps = strain_opt("max_shape_steps", 0);
    unsigned int substr_mat = strain_opt("substrate_material", 1);
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
    phys_reg_ID[0] = 1; // InP substrate
    phys_reg_ID[1] = 2; // AlInAs
    phys_reg_ID[2] = 3; // GaInAs
    phys_reg_ID[3] = 4; // AlInAs
    
    vector<unsigned int> BC_reg_ID(3);
    BC_reg_ID[0] = 1; // anode
    BC_reg_ID[1] = 2; // cathode

    
    unsigned int dim = 1;
    Mesh mesh(1);
    MeshTools::Generation::build_line(mesh, 
        200,
        0, 20,
        EDGE2);


    MeshData_elements meshdata(mesh);
    meshdata.enable_compatibility_mode();
    
    map<unsigned int, vector<unsigned int> > boundary_nodes;
    boundary_nodes[1] = vector<unsigned int>(1, 200);
    boundary_nodes[2] = vector<unsigned int>(1, 0);

    {
      map<const Elem*, vector<Number> > m;
      MeshBase::const_element_iterator it = mesh.elements_begin();
      const MeshBase::const_element_iterator end = mesh.elements_end();
      for ( ; it != end; ++it)
      {
        const Elem* elem = *it;

        double y = elem->centroid()(0);
        if (y < 5.0)
          m[elem] = vector<Number>(1, 1);
        else if (y < 10.0)
          m[elem] = vector<Number>(1, 2);
        else if (y < 15.0)
          m[elem] = vector<Number>(1, 3);
        else
          m[elem] = vector<Number>(1, 4);
      }
      meshdata.insert_elem_data(m);

    }

 
    mesh.print_info();


    /*****************************************************
     *
     * Setup of strain parameters
     *
     *****************************************************/

    bool periodicity[3];

    rotated_crystal cryst;
    // InAs - GaAs
    cryst.set_cryst_type("cub");
    cryst.set_xyz_mil_direction("x",  1,  0, 0);
    cryst.set_xyz_mil_direction("y",  0,  1, 0);
    cryst.set_xyz_mil_direction("z",  0,  0, 1);


    double x_al = al_content;
    double x_ga = ga_content;
    double x_in = 1 - al_content;
    
    Macrostrain::strain_param inp_strain;
    Macrostrain::strain_param inalas_strain;
    Macrostrain::strain_param gainas_strain;
    Piezoelectricity inp_piezo;
    Piezoelectricity inalas_piezo;
    Piezoelectricity gainas_piezo;

    std::map<unsigned int, Macrostrain::strain_param*> strain_params;
    strain_params[1] = &inp_strain;
    strain_params[2] = &inalas_strain;
    strain_params[3] = &gainas_strain;
    strain_params[4] = &inalas_strain;

    std::map<unsigned int, Piezoelectricity*> piezodata;
    piezodata[1] = &inp_piezo;
    piezodata[2] = &inalas_piezo;
    piezodata[3] = &gainas_piezo;
    piezodata[4] = &inalas_piezo;

    {

      inp_strain.C_tensor.set_moduli(101.1,  56.1,   45.6); // InP
      inalas_strain.C_tensor.set_moduli(
          alloy(83.290, 125.0, x_in),
          alloy(45.260, 53.4, x_in),
          alloy(39.590, 54.2, x_in)); //InAlAs
      gainas_strain.C_tensor.set_moduli(
          alloy(122.1, 83.290, x_ga),
          alloy(56.6, 45.260, x_ga),
          alloy(60.0, 39.590, x_ga)); //GaInAs

      inp_piezo.set_moduli(0.056); //InP
      inalas_piezo.set_moduli(alloy(-0.044, -0.015, x_in)); //InAs
      gainas_piezo.set_moduli(alloy(-0.16, -0.044, x_ga)); //GaAs  C/m^2
      inp_piezo.set_pyro_module(0.0);
      inalas_piezo.set_pyro_module(0.0);
      gainas_piezo.set_pyro_module(0.0);

      inp_strain.crystal = cryst;
      inalas_strain.crystal = cryst;
      gainas_strain.crystal = cryst;

      std::vector<int> x_dir(3);
      std::vector<int> y_dir(3);
      switch (growth_dir)
      {
        case 1:
          x_dir[0] = 1;  x_dir[1] = 0;  x_dir[2] = 0;
          y_dir[0] = 0;  y_dir[1] = 1;  y_dir[2] = 0;
          break;
        case 2:
          x_dir[0] = 0;  x_dir[1] = 1;  x_dir[2] = 1;
          y_dir[0] = 0;  y_dir[1] = 1;  y_dir[2] = -1;
          break;
        case 3:
          x_dir[0] = 1;  x_dir[1] = 1;  x_dir[2] = 1;
          y_dir[0] = 0;  y_dir[1] = 1;  y_dir[2] = -1;
          break;
      }

      inp_strain.crystal.set_lat_const(5.8697);
      inp_strain.crystal.calculate_lat_consts();
      inp_strain.crystal.calculate_rot_matrix(x_dir, y_dir);
      inp_strain.C_tensor.rotate_to_calc_system(
          inp_strain.crystal.RotMatrix);

      inalas_strain.crystal.set_lat_const(alloy(6.05830, 5.6611, x_in));
      inalas_strain.crystal.calculate_lat_consts();
      inalas_strain.crystal.calculate_rot_matrix(x_dir, y_dir);
      inalas_strain.C_tensor.rotate_to_calc_system(
          inalas_strain.crystal.RotMatrix);

      gainas_strain.crystal.set_lat_const(alloy(5.65325, 6.05830, x_ga));
      gainas_strain.crystal.calculate_lat_consts();
      gainas_strain.crystal.calculate_rot_matrix(x_dir, y_dir);
      gainas_strain.C_tensor.rotate_to_calc_system(
          gainas_strain.crystal.RotMatrix);

    }


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
    opt.mesh_input_file    = meshfile;

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

    double stress_value = strain_opt("stress", 0.0);
    map<unsigned int, double> stress_map;
    stress_map[1] = stress_value;


    Macrostrain strain_calculation(opt, mesh);

    strain_calculation.define_substrate_bc(2);
    strain_calculation.define_BC_map(boundary_nodes);
    strain_calculation.define_stress_value(stress_map);

    strain_calculation.assign_mesh_data(meshdata);
    strain_calculation.define_strain_parameters(strain_params);
    strain_calculation.define_piezo_moduli(piezodata);
    //strain_calculation.define_external_stress(stress_vector_in,
    //    boundary_nodes);


    /*****************************************************/


    /****************************************************************
     *
     * Drift diffusion definitions
     *
     ****************************************************************/

    StrainedSemiconductorModel inp(&strain_calculation);    
    inp.set_data_file("materials/InP.dat");

    if (statistics == "FD")
      inp.set_statistics(TiberCad::FERMIDIRAC);
    else
      inp.set_statistics(TiberCad::BOLTZMANN);
    
    if (!solve_strain)
      inp.ignore_strain();

    // InP
    //inp.add_recombination_model(SRH);
    //inp.set_SRH_parameters(1e-7, 1e-7);
    //inp.add_recombination_model(DIRECT);
    //inp.set_direct_rec_parameters(1e-15);

    // AlInAs
    StrainedSemiconductorModel alinas(inp);
    alinas.set_data_file("materials/InAs.dat");
    
    // InGaAs
    StrainedSemiconductorModel gainas(inp);
    gainas.set_data_file("materials/InAs.dat");

    //gainas.set_SRH_parameters(1e-7, 1e-7);
    //gainas.add_recombination_model(DIRECT);
    //gainas.set_direct_rec_parameters(1e-15);

    //inp.set_n_dopant(Dopant(n_doping, 0.025, 2));
    //alinas.set_p_dopant(Dopant(p_doping, 0.01, 4));

    Dummy d;
    inp.read_database(d);
    alinas.read_database(d);
    gainas.read_database(d);

    alinas.build_alloy("materials/AlAs.dat",
        "materials/AlInAs_bow.dat", x_al);
    gainas.build_alloy("materials/GaAs.dat",
        "materials/GaInAs_bow.dat", x_ga);

    ElementData element_data;
    {
      MeshData::const_elem_data_iterator it = meshdata.elem_data_begin();
      const MeshData::const_elem_data_iterator end =
        meshdata.elem_data_end();
      for ( ; it != end; ++it)
      {
        const Elem* elem = it->first;

        // every element needs to have a material assigned
        assert(meshdata.has_data(elem));

        int id = (int) meshdata(elem);

        switch (id)
        {
          case 1:
            element_data.set_data(elem, &inp);
            break;
          case 2:
            element_data.set_data(elem, &alinas);
            break;
          case 3:
            element_data.set_data(elem, &gainas);
            break;
          case 4:
            element_data.set_data(elem, &alinas);
            break;
        }
      }
    }



    // solve strain
    if (solve_strain)
    {
      cout << "Solving strain... \n" << flush;
      strain_calculation.solve();
      strain_calculation.output_piezo("Piezo.gmv");
    }


    if (solve_strain)
    {
      inp.include_strain();
      alinas.include_strain();
      gainas.include_strain();
    }


    OhmicContact anode("anode");
    OhmicContact cathode("cathode");

    BoundaryData boundary_data;
    {
      map<unsigned int, vector<unsigned int> >::const_iterator it =
        boundary_nodes.begin();
      const map<unsigned int, vector<unsigned int> >::const_iterator end =
        boundary_nodes.end();

      for ( ; it != end; ++it)
      {
        const vector<unsigned int>& nodes = it->second;

        switch (it->first)
        {
          case 1:
            //set_boundary(boundary_data, nodes, &anode, mesh);
            break;
          case 2:
            //set_boundary(boundary_data, nodes, &cathode, mesh);
            break;
        }
      }
    }

    DD::Device device(&mesh, &element_data, &boundary_data);
    bool device_integrity = device.check_integrity();
    if (device_integrity)
      cout << "Device ok.\n\n";
    else
    {
      cout << "Device bad.\n\n";
      return 1;
    }

    mesh.print_info();
    
  
    DriftDiffusion dd(&device);

    DriftDiffusion::Options& params = dd.get_options();
    params.max_refinement_steps = dd_max_r_steps;
    params.solver_params.nonlinear_max_iterations = 25;
    params.solver_params.linear_max_iterations = dd_lin_max_it;
    params.solver_params.nonlinear_tolerance = dd_nonlin_rtol;
    params.solver_params.nonlinear_abs_tolerance = dd_nonlin_atol;
    params.solver_params.ls_maxstep = 0.05;
    params.solver_params.linear_tolerance = dd_lin_rtol;
    params.solver_params.linear_abs_tolerance = dd_lin_atol;
    params.integration_order =
      static_cast<libMeshEnums::Order>(integration_order);

    params.solver_params.pc_type = PCILU;
    params.local_scaling = true;
    

    if (method == "GUMMEL")
    {
      params.solver_method = DriftDiffusion::GUMMEL;
      params.max_gummel_iterations = 2;
    }
    else
      params.solver_method = DriftDiffusion::NEWTON;

    if (fully)
    {
      cout << "Solving for electrons and holes.\n" << flush;
      params.coupling = FULLYCOUPLED;
    }
    else
    {
      cout << "Solving for electrons only.\n" << flush;
      params.coupling = POISSON | ELECTRONS;
    }


    if (approx_order == "FIRST")
      params.approximation_order = FIRST;
    else if (approx_order == "SECOND")
      params.approximation_order = SECOND;

    // mesh drawn in um
    params.mesh_units = mesh_units;

    dd.enable_mesh_refinement();

    /*****************************************************/

    


    cout << "Solving drift-diffusion... \n" << flush;

    dd.set_simulation_voltage("cathode", 0.0);
    dd.set_simulation_voltage("anode", 0.0);
    
    cout << "Solving equilibrium... " << flush;
    dd.solve();
    const Scaling& sc = dd.get_scaling();
    cout << "Scaling parameters:\n";
    cout << "     phi0: " << sc.get_potential_scaling() << "\n";
    cout << "     x0  : " << sc.get_length_scaling() << "\n";
    cout << "     mu0 : " << sc.get_mobility_scaling() << "\n";
    cout << "     C0  : " << sc.get_density_scaling() << "\n";
    cout << "     t0  : " << sc.get_time_scaling() << "\n";
    cout << "     R0  : " << sc.get_recombination_scaling() << "\n\n";

    cout << "InP:\n";
    inp.print_info();
    cout << "AlInAs:\n";
    alinas.print_info();
    cout << "GaInAs:\n";
    gainas.print_info();

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

    params.solver_params.nonlinear_max_iterations = dd_nonlin_max_it;
    params.solver_params.ls_maxstep = dd_nonlin_ls_maxstep;

    cout << "\nBegin sweep...\n" << flush;
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
    for (int i = 0; i < voltage_steps; i++)
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

      const map<const ElectricalContact*, double>& curr =
        dd.get_boundary_currents();
      file << *it << "  "
           << (*curr.find(&cathode)).second << "  "
           << (*curr.find(&anode)).second << "  "
           << dd.get_artificial_boundary_current() << "\n" << flush;
      cerr << "    I = " << (*curr.find(&cathode)).second << " A/cm^2\n";
    }


    it = zero;
    if (it != voltages.begin())
    {
      bool restart = true;
      do
      {
        --it;
        dd.set_simulation_voltage("anode", *it);
        cout << " Solving U = " << *it << " V ... " << flush;
        dd.solve(restart);
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

        restart = false;
        const map<const ElectricalContact*, double>& curr =
          dd.get_boundary_currents();
        file << *it << "  "
          << (*curr.find(&cathode)).second << "  "
          << (*curr.find(&anode)).second << "  "
          << dd.get_artificial_boundary_current() << "\n" << flush;
        cerr << "    I = " << (*curr.find(&cathode)).second << " A/cm^2\n";
      }
      while (it != voltages.begin());
    }
    file.close();
  }

  return libMesh::close();
}


void set_boundary(BoundaryData& data, const vector<unsigned int>& nodes,
    ElectricalContact* desc, const Mesh& mesh)
{
  vector<int>::const_iterator n_it;
  const vector<unsigned int>::const_iterator n_begin = nodes.begin();
  const vector<unsigned int>::const_iterator n_end = nodes.end();

  Mesh::const_element_iterator el = mesh.elements_begin();
  const Mesh::const_element_iterator el_end = mesh.elements_end();
  for ( ; el != el_end; ++el)
  {
    Elem* elem = *el;

    int n_sides = elem->n_sides();
    for (int s = 0; s < n_sides; s++)
    {
      if (elem->neighbor(s) == NULL)
      {
        if (find(n_begin, n_end, elem->node(s)) != n_end)
          data.set_data(BoundaryData::ElementSide(elem, s), desc);
      }
    }
  }
}

