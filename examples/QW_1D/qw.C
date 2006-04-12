// the following _HAS_ to be included first
#include "Read_MSH.h"

#include "ElementData.h"
#include "BoundaryDescriptor.h"
#include "BoundaryData.h"
#include "DDevice.h"
#include "DriftDiffusion.h"
#include "SemiconductorModel.h"

#include "macrostrain.h"
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

void setup_boundary_desc(BoundaryDescriptor& desc,
    DriftDiffusionProperties& sc_model);

void set_boundary(BoundaryData& data, const vector<unsigned int>& nodes,
    BoundaryDescriptor& desc, const Mesh& mesh);

void set_boundary(BoundaryData& data, double x,
    BoundaryDescriptor& desc, const Mesh& mesh);


double alloy(double a, double b, double xa, double bowing = 0.0)
{
  return b + (a - b) * xa - bowing * xa *  (1 - xa);
};


double schottky_barrier;


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

    int mat_sys = input_file("material_system",0);
    double in_content = input_file("in_content", 0.15);
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
        103,
        -51.5, 51.5,
        EDGE2);

    MeshData_elements meshdata(mesh);
    meshdata.enable_compatibility_mode();
    
    map<unsigned int, vector<unsigned int> > boundary_nodes;
    boundary_nodes[1] = vector<unsigned int>(1, 104);
    boundary_nodes[2] = vector<unsigned int>(1, 1);

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

 
    mesh.print_info();


    /*****************************************************
     *
     * Setup of strain parameters
     *
     *****************************************************/
    BC_region_type substrate_region;

    std::vector<BC_region_type> material_region;

    bool periodicity[3];

    stiffness C1;
    rotated_crystal cryst;

    // AlGaN - GaN
    cryst.set_cryst_type("hex");
    cryst.set_xyz_mil_direction("y", 1,  0, -1, 0) ;
    cryst.set_xyz_mil_direction("z", 1,  -2, 1, 0) ;
    cryst.set_xyz_mil_direction("x",  0,  0, 0,  1) ;

    Piezoelectricity piezo1;

    std::vector<stiffness> C_tensor(number_of_regions, C1);
    std::vector<rotated_crystal> crystal(number_of_regions, cryst);
    std::vector<Piezoelectricity> piezo_data(number_of_regions, piezo1);

    std::vector<int> x_dir(4);
    std::vector<int> y_dir(4);

    y_dir[0] =  1;  y_dir[1] = 0;  y_dir[2] = -1;  y_dir[3] =  0;
    x_dir[0] =  0;  x_dir[1] = 0;  x_dir[2] =  0;  x_dir[3] =  1;


    // GaN
    C_tensor[0].set_moduli(390.0, 145.0, 106.0, 398.0, 105.0);
    C_tensor[2].set_moduli(390.0, 145.0, 106.0, 398.0, 105.0);
    //
    piezo_data[0].set_moduli(1.27, -0.35, -0.3); // C/m^2
    piezo_data[2].set_moduli(1.27, -0.35, -0.3); // C/m^2
    //
    piezo_data[0].set_pyro_module(-0.034);
    piezo_data[2].set_pyro_module(-0.034);
    //
    crystal[0].set_lat_const(0.3189, 0.5185);
    crystal[2].set_lat_const(0.3189, 0.5185);
    crystal[0].calculate_lat_consts();
    crystal[2].calculate_lat_consts();
    crystal[0].calculate_rot_matrix(x_dir, y_dir); 
    crystal[2].calculate_rot_matrix(x_dir, y_dir); 
    C_tensor[0].rotate_to_calc_system(crystal[0].RotMatrix);
    C_tensor[2].rotate_to_calc_system(crystal[2].RotMatrix);
    
    // InGaN
    // first value: InN, second value: GaN
    C_tensor[1].set_moduli(
        alloy(223.0, 390.0, in_content),
        alloy(115.0, 145.0, in_content),
        alloy(92.0, 106.0, in_content),
        alloy(224.0, 398.0, in_content),
        alloy(48.0, 105.0, in_content));

    piezo_data[1].set_moduli(
        alloy(0.97, 1.27, in_content),
        alloy(-0.57, -0.35, in_content),
        alloy(-0.4, -0.3, in_content));
    piezo_data[1].set_pyro_module(
        alloy(-0.042, -0.034, in_content, -0.037));

    crystal[1].set_lat_const(
        alloy(0.3545, 0.3189, in_content),
        alloy(0.5703, 0.5185, in_content));
    crystal[1].calculate_lat_consts();
    crystal[1].calculate_rot_matrix(x_dir, y_dir);
    C_tensor[1].rotate_to_calc_system(crystal[1].RotMatrix);

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

    for(int i = 0; i<=2; i++)
      substrate_region.coord_min[i] = min_coord_substrate[i];
    for(int i = 0; i<=2; i++)
      substrate_region.coord_max[i] = max_coord_substrate[i];


    vector<double> fp(3);

    fp[0] = fp[1] = fp[2] = 0.0;
    opt.fixed_point1 = fp;
    fp[0] = 100;
    opt.fixed_point2 = fp;
    fp[1] = 0.0; fp[0] = fp[2] = 1.0;
    opt.fixed_point3 = fp;

    double stress_value = strain_opt( "stress",0.0);
    vector<external_stress>  stress_vector_in(1);
    stress_vector_in[0].bc_region_number = 1;
    stress_vector_in[0].stress_value = stress_value;


    Macrostrain strain_calculation(opt, mesh);

    strain_calculation.assign_mesh_data(meshdata);
    strain_calculation.define_strain_parameters(C_tensor, crystal);
    strain_calculation.define_substrate_region(substrate_region);
    strain_calculation.define_piezo_moduli(piezo_data);
    //strain_calculation.define_external_stress(stress_vector_in,
    //    boundary_nodes);


    /*****************************************************/



    /****************************************************************
     *
     * Drift diffusion definitions
     *
     ****************************************************************/

    StrainedSemiconductorModel n_gan(&strain_calculation);    
    if (!solve_strain)
      n_gan.ignore_strain();

    if (statistics == "FD")
      n_gan.set_statistics(TiberCad::FERMIDIRAC);
    else
      n_gan.set_statistics(TiberCad::BOLTZMANN);

    // n-GaN
    n_gan.set_relative_permittivity(9.5);
    n_gan.set_valence_band_properties(-0.72, 1.6, 200);
    n_gan.set_conduction_band_properties(2.789, 0.23, 1000);
    n_gan.add_recombination_model(SRH);
    n_gan.set_SRH_parameters(1e-8, 1e-8);
    n_gan.add_recombination_model(DIRECT);
    n_gan.set_direct_rec_parameters(1e-12);

    // p-GaN
    StrainedSemiconductorModel p_gan(n_gan);
    
    // InGaN
    StrainedSemiconductorModel ingan(n_gan);
    ingan.set_relative_permittivity(alloy(14.0, 9.5, in_content));
    ingan.set_valence_band_properties(
        alloy(-0.462, -0.72, in_content),
        alloy(1.65, 1.6, in_content),
        alloy(200, 200, in_content));
    ingan.set_conduction_band_properties(
        alloy(0.318, 2.789, in_content, 1.4),
        alloy(0.11, 0.23, in_content),
        alloy(1000, 1000, in_content));
    ingan.set_SRH_parameters(1e-8, 1e-8);
    ingan.add_recombination_model(DIRECT);
    ingan.set_direct_rec_parameters(1e-12);
    //ingan.set_relative_permittivity(9.5);
    //ingan.set_valence_band_properties(-0.72, 1.6, 200);
    //ingan.set_conduction_band_properties(2.789, 0.23, 1000);

    n_gan.set_n_dopant(Dopant(n_doping, 0.025, 2));
    p_gan.set_p_dopant(Dopant(p_doping, 0.01, 4));
    //ingan.set_p_dopant(Dopant(1e17, 0.01, 4));


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
            element_data.set_data(elem, &n_gan);
            break;
          case 2:
            element_data.set_data(elem, &ingan);
            break;
          case 3:
            element_data.set_data(elem, &p_gan);
            break;
        }
      }
    }


///*
    {
      MeshRefinement ref(mesh);
      ref.uniformly_refine();
      ref.uniformly_refine();
      ref.uniformly_refine();
      ref.uniformly_refine();
      //ref.uniformly_refine();

      //for (int i = 0; i < 3; i++)
      for (int i = 0; i < 0; i++)
      {
        {
          MeshBase::element_iterator it = mesh.elements_begin();
          const MeshBase::element_iterator end = mesh.elements_end();
          for ( ; it != end; ++it)
          {
            Elem* elem = *it;

            double y = elem->centroid()(0);
            if ((y > -2) && (y < 2))
              elem->set_refinement_flag(Elem::REFINE);
          }
          ref.refine_elements();
        }
        /*
        {
          MeshBase::element_iterator it = mesh.elements_begin();
          const MeshBase::element_iterator end = mesh.elements_end();
          for ( ; it != end; ++it)
          {
            Elem* elem = *it;

            double y = elem->centroid()(0);
            if ((y > -10.0) && (y < 10.0))
              elem->set_refinement_flag(Elem::REFINE);
          }
          ref.refine_elements();
        }
        {
          MeshBase::element_iterator it = mesh.elements_begin();
          const MeshBase::element_iterator end = mesh.elements_end();
          for ( ; it != end; ++it)
          {
            Elem* elem = *it;

            double y = elem->centroid()(0);
            if ((y > -2.5) && (y < 2.5))
              elem->set_refinement_flag(Elem::REFINE);
          }
          ref.refine_elements();
        }
        */
      }
    }
//*/

    // solve strain
    if (solve_strain)
    {
      cout << "Solving strain... \n" << flush;
      strain_calculation.solve();
      strain_calculation.output_piezo("Piezo.gmv");
    }


    const Elem* elem = meshdata.elem_data_begin()->first;
    n_gan.reinit(elem);
    p_gan.reinit(elem);
    ingan.reinit(elem);

    cout << "  n-GaN:" << endl;
    n_gan.calculate_equilibrium_properties();
    cout << "       " << n_gan.get_equilibrium_fermi_level() << " eV, "
      << "ni = " << n_gan.get_intrinsic_density()
      << " n0 = " << n_gan.get_equilibrium_electron_density()
      << " p0 = " << n_gan.get_equilibrium_hole_density() << endl;

    cout << "  p-GaN:" << endl;
    p_gan.calculate_equilibrium_properties();
    cout << "       " << p_gan.get_equilibrium_fermi_level() << " eV, "
      << "ni = " << p_gan.get_intrinsic_density()
      << " n0 = " << p_gan.get_equilibrium_electron_density()
      << " p0 = " << p_gan.get_equilibrium_hole_density() << endl;

    cout << "  InGaN:" << endl;
    ingan.calculate_equilibrium_properties();
    cout << "       " << ingan.get_equilibrium_fermi_level() << " eV, "
      << "ni = " << ingan.get_intrinsic_density()
      << " n0 = " << ingan.get_equilibrium_electron_density()
      << " p0 = " << ingan.get_equilibrium_hole_density() << endl;




    BoundaryDescriptor anode("anode");
    BoundaryDescriptor cathode("cathode");
    setup_boundary_desc(anode, p_gan);
    setup_boundary_desc(cathode, n_gan);


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
            set_boundary(boundary_data, nodes, anode, mesh);
            break;
          case 2:
            set_boundary(boundary_data, nodes, cathode, mesh);
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
    params.solver_params.ls_maxstep = 0.1;
    params.solver_params.linear_tolerance = dd_lin_rtol;
    params.solver_params.linear_abs_tolerance = dd_lin_atol;
    params.integration_order =
      static_cast<libMeshEnums::Order>(integration_order);

    params.solver_params.pc_type = PCILU;
    params.artificial_drift = true;
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

      const map<const BoundaryDescriptor*, double>& curr =
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
        const map<const BoundaryDescriptor*, double>& curr =
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



void setup_boundary_desc(BoundaryDescriptor& desc,
      DriftDiffusionProperties& sc)
{

  std::vector<double> coeff(3, 0);

  coeff[0] = 1.0;
  if (desc.get_id() == "anode")
  {
    desc.set_coefficients("fermi_h", coeff);
    coeff[0] = 0.0; coeff[1] = 1.0;
    desc.set_coefficients("fermi_e", coeff);
  }
  else
  {
    desc.set_coefficients("fermi_e", coeff);
    coeff[0] = 0.0; coeff[1] = 1.0;
    desc.set_coefficients("fermi_h", coeff);
  }

  coeff[1] = 0.0; coeff[0] = 1.0;
  coeff[2] = sc.get_equilibrium_fermi_level();
  desc.set_coefficients("potential", coeff);
}


void set_boundary(BoundaryData& data, const vector<unsigned int>& nodes,
    BoundaryDescriptor& desc, const Mesh& mesh)
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
        if (find(n_begin, n_end, elem->node(s) + 1) != n_end)
          data.set_data(BoundaryData::ElementSide(elem, s), &desc);
      }
    }
  }
}

void set_boundary(BoundaryData& data, double x,
    BoundaryDescriptor& desc, const Mesh& mesh)
{

  MeshBase::const_element_iterator el = mesh.active_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_elements_end();
  while (el != end_el)
  {
    const Elem* elem = *el;

    for (unsigned int n = 0; n < elem->n_nodes(); n++)
    {
      if (elem->neighbor(n) == NULL) {
        double xcoord = (*elem->get_node(n))(0);
        
        if (std::fabs(xcoord - x) < 1e-15)
        {
          cerr << "Boundary: " << xcoord << " " << elem->node(n) << "\n";
          data.set_data(BoundaryData::ElementSide(elem, n), &desc);
        }
      }
    }

    ++el;
  }
}

