// the following _HAS_ to be included first
#include "Read_MSH.h"

#include "ElementData.h"
#include "OhmicContact.h"
#include "SchottkyContact.h"
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
#include "gmv_io.h"

#include <algorithm>

using namespace std;
using namespace DriftDiffusionDefs;

void set_boundary(BoundaryData& data, const vector<unsigned int>& nodes,
    ElectricalContact* desc, const Mesh& mesh);

void sweep_drain(double stop, int steps, DriftDiffusion& dd, double vg,
    bool restart = false);

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
    const double mesh_units = input_file("mesh_units", 1e-4);

    int fully = input_file("fully_coupled", 1);

    const string approx_order =
      input_file("approximation_order", "FIRST");
    unsigned int dim = input_file("dimension",3);            
    int mat_sys = input_file("material_system",0);
    double al_content = input_file("al_content", 1.0);
    int growth_dir = input_file("growth_dir",3);

    // drift-diffusion options
    GetPot dd_opt("dd.in");
    int characteristic = dd_opt("characteristic", 0);
    double n_doping = dd_opt("n_doping", 0.0);

    double vds_stop = dd_opt("vds_stop", 0.0);
    unsigned int vds_steps = dd_opt("vds_steps", 1);
    double vg_start = dd_opt("vg_start", 0.0);
    double vg_stop = dd_opt("vg_stop", 0.0);
    unsigned int vg_steps = dd_opt("vg_steps", 1);

    double vg = dd_opt("vg", 0.0);
    double vds = dd_opt("vds", 0.0);

    double schottky_barrier = dd_opt("schottky_barrier", 1.2);

    const string method = dd_opt("simulation_method", "NEWTON");
    const string statistics = dd_opt("statistics", "B");
    double dd_nonlin_rtol = dd_opt("nonlinear_tolerance", 1e-9);
    double dd_nonlin_atol = dd_opt("nonlinear_abs_tolerance", 1e-18);
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

    
    unsigned int number_of_regions = 4;
    vector<unsigned int> phys_reg_ID(number_of_regions);
    phys_reg_ID[0] = 1; // AlN
    phys_reg_ID[1] = 2; // GaN (channel)
    phys_reg_ID[2] = 3; // AlN doped
    phys_reg_ID[3] = 4; // GaN doped
    
    vector<unsigned int> BC_reg_ID(4);
    BC_reg_ID[0] = 1; // gate and stress
    BC_reg_ID[1] = 2; // source
    BC_reg_ID[2] = 3; // drain
    BC_reg_ID[3] = 4; // substrate

    
    Mesh mesh(dim);
    MeshData_elements meshdata(mesh);
    meshdata.enable_compatibility_mode();
    
    Read_MSH readmesh(meshfile, phys_reg_ID, BC_reg_ID, dim, mesh, meshdata);
    map<unsigned int, vector<unsigned int> > boundary_nodes;
    readmesh.get_BC_data(boundary_nodes);
    
    mesh.print_info();


    EquationSystems eqsys(mesh);

    
    /*****************************************************
     *
     * Setup of strain parameters
     *
     *****************************************************/
    bool periodicity[3];

    stiffness C1;
    rotated_crystal cryst;

    if (mat_sys == 1)
    {
      // InAs - GaAs
      cryst.set_cryst_type("cub");
      cryst.set_xyz_mil_direction("x",  1,  0, 0);
      cryst.set_xyz_mil_direction("y",  0,  1, 0);
      cryst.set_xyz_mil_direction("z",  0,  0, 1);
    }
    else
    {
      // AlGaN - GaN
      cryst.set_cryst_type("hex");
      cryst.set_xyz_mil_direction("z", 1,  0, -1, 0) ;
      cryst.set_xyz_mil_direction("x", 1,  -2, 1, 0) ;
      cryst.set_xyz_mil_direction("y",  0,  0, 0,  1) ;
    }

    Piezoelectricity piezo1;

    std::vector<stiffness> C_tensor(number_of_regions, C1);
    std::vector<rotated_crystal> crystal(number_of_regions, cryst);
    std::vector<Piezoelectricity> piezo_data(number_of_regions, piezo1);


    if (mat_sys == 1)
    {
      // GaAs - InAs
      cout << "Using GaAs - InAs\n" << flush;

      C_tensor[0].set_moduli(122.1,  56.6,   60.0); //GaAs
      C_tensor[1].set_moduli(83.290,  45.260,  39.590); //InAs
      C_tensor[2].set_moduli(122.1,  56.6,   60.0); //GaAs doped
      C_tensor[3].set_moduli(83.290,  45.260,  39.590); //InAs doped

      piezo_data[0].set_moduli(0.16); //GaAs  C/m^2
      piezo_data[1].set_moduli(0.044); //InAs
      piezo_data[2].set_moduli(0.16); //GaAs doped
      piezo_data[3].set_moduli(0.044); //InAs doped
      piezo_data[0].set_pyro_module(0.0);
      piezo_data[1].set_pyro_module(0.0);
      piezo_data[2].set_pyro_module(0.0);
      piezo_data[3].set_pyro_module(0.0);


      std::vector<int> x_dir(3);
      std::vector<int> y_dir(3);
      // growth direction
      y_dir[0] = growth_dir;  y_dir[1] = 1;  y_dir[2] =   1;
      x_dir[0] = 0;  x_dir[1] = 1;  x_dir[2] =  -1;

      crystal[0].set_lat_const(5.65325);
      crystal[0].calculate_lat_consts();
      crystal[0].calculate_rot_matrix(x_dir, y_dir);
      C_tensor[0].rotate_to_calc_system(crystal[0].RotMatrix);

      crystal[1].set_lat_const(6.05830);
      crystal[1].calculate_lat_consts();
      crystal[1].calculate_rot_matrix(x_dir, y_dir);
      C_tensor[1].rotate_to_calc_system(crystal[1].RotMatrix);

      crystal[2].set_lat_const(5.65325);
      crystal[2].calculate_lat_consts();
      crystal[2].calculate_rot_matrix(x_dir, y_dir);
      C_tensor[2].rotate_to_calc_system(crystal[2].RotMatrix);

      crystal[3].set_lat_const(6.05830);
      crystal[3].calculate_lat_consts();
      crystal[3].calculate_rot_matrix(x_dir, y_dir);
      C_tensor[3].rotate_to_calc_system(crystal[3].RotMatrix);
    }
    else
    {
      // AlGaN-GaN
      cout << "Using AlGaN - GaN (with "
        << al_content << " Al content)\n" << flush;
      
      C_tensor[1].set_moduli(390.0, 145.0, 106.0, 398.0, 105.0); //GaN
      C_tensor[0].set_moduli(
          alloy(396.0, 390.0, al_content),
          alloy(137.0, 145.0, al_content),
          alloy(108.0, 106.0, al_content),
          alloy(373.0, 398.0, al_content),
          alloy(116.0, 105.0, al_content)); //AlN

      piezo_data[1].set_moduli(1.27, -0.35, -0.3); //GaN  C/m^2
      piezo_data[1].set_pyro_module(-0.034); //GaN
      piezo_data[0].set_moduli(
          alloy(1.79, 1.27, al_content),
          alloy(-0.5, -0.35, al_content),
          alloy(-0.48, -0.3, al_content)); //AlN
      piezo_data[0].set_pyro_module(
          alloy(-0.09, -0.034, al_content, -0.021)); //AlN



      std::vector<int> x_dir(4);
      std::vector<int> y_dir(4);

      //x_dir[0] = 1;  x_dir[1] = 0;  x_dir[2] =  -1;  x_dir[3] = 0;
      x_dir[0] = 1;  x_dir[1] = -2;  x_dir[2] = 1;  x_dir[3] = 0;
      y_dir[0] = 0;  y_dir[1] = 0;  y_dir[2] =  0;  y_dir[3] = 1;

      crystal[1].set_lat_const(0.3189, 0.5185); //GaN
      crystal[1].calculate_lat_consts();
      crystal[1].calculate_rot_matrix(x_dir, y_dir); 
      C_tensor[1].rotate_to_calc_system(crystal[1].RotMatrix);

      crystal[0].set_lat_const(
          alloy(0.3112, 0.3189, al_content),
          alloy(0.4982, 0.5185, al_content)); //AlN
      crystal[0].calculate_lat_consts();
      crystal[0].calculate_rot_matrix(x_dir, y_dir);
      C_tensor[0].rotate_to_calc_system(crystal[0].RotMatrix);

      // doped parts
      C_tensor[3].set_moduli(390.0, 145.0, 106.0, 398.0, 105.0); //GaN doped
      C_tensor[2].set_moduli(
          alloy(396.0, 390.0, al_content),
          alloy(137.0, 145.0, al_content),
          alloy(108.0, 106.0, al_content),
          alloy(373.0, 398.0, al_content),
          alloy(116.0, 105.0, al_content)); //AlN doped

      piezo_data[3].set_moduli(1.27, -0.35, -0.3); //GaN  C/m^2
      piezo_data[3].set_pyro_module(-0.034); //GaN
      piezo_data[2].set_moduli(
          alloy(1.79, 1.27, al_content),
          alloy(-0.5, -0.35, al_content),
          alloy(-0.48, -0.3, al_content)); //AlN
      piezo_data[2].set_pyro_module(
          alloy(-0.09, -0.034, al_content, -0.021)); //AlN

      crystal[3].set_lat_const(0.3189, 0.5185); //GaN
      crystal[3].calculate_lat_consts();
      crystal[3].calculate_rot_matrix(x_dir, y_dir); 
      C_tensor[3].rotate_to_calc_system(crystal[3].RotMatrix);

      crystal[2].set_lat_const(
          alloy(0.3112, 0.3189, al_content),
          alloy(0.4982, 0.5185, al_content)); //AlN
      crystal[2].calculate_lat_consts();
      crystal[2].calculate_rot_matrix(x_dir, y_dir);
      C_tensor[2].rotate_to_calc_system(crystal[2].RotMatrix);
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
    fp[1] = -0.5;
    opt.fixed_point2 = fp;
    fp[1] = 0.0; fp[0] = fp[2] = 1.0;
    opt.fixed_point3 = fp;

    double stress_value = strain_opt( "stress",0.0);
    map<unsigned int, double> stress_map;
    stress_map[1] = stress_value;

   
    Macrostrain::strain_param GaN_strain;
    Macrostrain::strain_param AlGaN_strain;
    Macrostrain::strain_param GaN_doped_strain;
    Macrostrain::strain_param AlGaN_doped_strain;


    GaN_strain.crystal = crystal[1] ;
    AlGaN_strain.crystal = crystal[0] ;
    GaN_doped_strain.crystal = crystal[3];
    AlGaN_doped_strain.crystal = crystal[2];

    GaN_strain.C_tensor = C_tensor[1];
    AlGaN_strain.C_tensor = C_tensor[0];
    GaN_doped_strain.C_tensor = C_tensor[3];
    AlGaN_doped_strain.C_tensor = C_tensor[2];
    

    std::map<unsigned int, Macrostrain::strain_param*> strain_params;

    strain_params[1] = &GaN_strain;
    strain_params[2] = &AlGaN_strain;
    strain_params[3] = &GaN_strain;
    strain_params[4] = &AlGaN_strain;
    

    std::map<unsigned int, Piezoelectricity*> piezodata;
    piezodata[1] = &piezo_data[1];
    piezodata[2] = &piezo_data[0];
    piezodata[3] = &piezo_data[3];
    piezodata[4] = &piezo_data[2];



    Macrostrain strain_calculation(opt, eqsys, "strainsys");

    strain_calculation.define_substrate_bc(4);
    strain_calculation.define_BC_map(boundary_nodes);
    strain_calculation.define_stress_value(stress_map);

    strain_calculation.assign_mesh_data(meshdata);

    strain_calculation.define_strain_parameters(strain_params);

    strain_calculation.define_piezo_moduli(piezodata);


    /*****************************************************/



    /****************************************************************
     *
     * Drift diffusion definitions
     *
     ****************************************************************/
    StrainedSemiconductorModel barrier(&strain_calculation);    

    if (statistics == "FD")
      barrier.set_statistics(TiberCad::FERMIDIRAC);
    else
      barrier.set_statistics(TiberCad::BOLTZMANN);
    

  
    double bg_doping = 1e16;
   
    barrier.add_dopant(new Dopant(bg_doping, 0.025, 2, Dopant::N_TYPE));


    StrainedSemiconductorModel channel(barrier);
    channel.set_n_dopant(Dopant(0, 0.025, 2));

    Dummy d;
    if (mat_sys == 1)
    {
      // GaAs
      barrier.set_SRH_parameters(1e-8, 1e-8);

      // InAs
      channel.set_SRH_parameters(1e-8, 1e-8);
    }
    else
    {
      // AlN
      barrier.set_data_file("materials/GaN.dat");
      barrier.set_SRH_parameters(1e-8, 1e-8);
      barrier.set_mobilities(1000, 300);
      barrier.read_database(d);
      barrier.build_alloy("materials/AlN.dat",
          "materials/AlGaN_bow.dat", al_content);

      // GaN
      channel.set_data_file("materials/GaN.dat");
      channel.set_SRH_parameters(1e-8, 1e-8);
      channel.set_mobilities(1000, 300);
      channel.read_database(d);
    }

    StrainedSemiconductorModel channel_doped(channel);
    channel_doped.set_data_file("materials/GaN.dat");
    channel_doped.read_database(d);
    StrainedSemiconductorModel barrier_doped(barrier);
    barrier_doped.set_data_file("materials/GaN.dat");
    barrier_doped.read_database(d);
    barrier_doped.build_alloy("materials/AlN.dat",
        "materials/AlGaN_bow.dat", al_content);
    barrier_doped.set_n_dopant(Dopant(n_doping, 0.025, 2));
    channel_doped.set_n_dopant(Dopant(n_doping, 0.025, 2));
    barrier_doped.set_SRH_parameters(1e-10, 1e-10);
    channel_doped.set_SRH_parameters(1e-10, 1e-10);


    if (fully)
    {
      barrier.set_coupling_type(BOTH);
      channel.set_coupling_type(BOTH);
      barrier_doped.set_coupling_type(BOTH);
      channel_doped.set_coupling_type(BOTH);
    }
    else
    {
      barrier.set_coupling_type(ELECTRONS);
      channel.set_coupling_type(ELECTRONS);
      barrier_doped.set_coupling_type(ELECTRONS);
      channel_doped.set_coupling_type(ELECTRONS);
    }

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
            element_data.set_data(elem, &barrier);
            break;
          case 2:
            element_data.set_data(elem, &channel);
            break;
          case 3:
            element_data.set_data(elem, &barrier_doped);
            break;
          case 4:
            element_data.set_data(elem, &channel_doped);
            break;
        }
      }
    }

    // solve strain
    {
      cout << "Solving strain... \n" << flush;
      strain_calculation.solve();
      ostringstream f;
      f.precision(5);
      f << "Piezo.gmv";
      strain_calculation.output_piezo(f.str());
    }


    OhmicContact source("source");
    SchottkyContact gate("gate");
    gate.set_schottky_barrier(schottky_barrier);
    OhmicContact drain("drain");

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
            set_boundary(boundary_data, nodes, &gate, mesh);
            break;
          case 2:
            set_boundary(boundary_data, nodes, &source, mesh);
            break;
          case 3:
            set_boundary(boundary_data, nodes, &drain, mesh);
            break;
          case 4:
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
    
  
    DriftDiffusion dd(&device);

    DriftDiffusion::Options& params = dd.get_options();
    params.max_refinement_steps = dd_max_r_steps;
    params.solver_params.nonlinear_max_iterations = dd_nonlin_max_it;
    params.solver_params.linear_max_iterations = dd_lin_max_it;
    params.solver_params.nonlinear_tolerance = dd_nonlin_rtol;
    params.solver_params.nonlinear_abs_tolerance = dd_nonlin_atol;
    params.solver_params.ls_maxstep = 0.2;
    params.solver_params.linear_tolerance = dd_lin_rtol;
    params.solver_params.linear_abs_tolerance = dd_lin_atol;
    params.integration_order =
      static_cast<libMeshEnums::Order>(integration_order);

    params.solver_params.pc_type = PCILU;
    params.local_scaling = true;
    //params.artificial_drift = false;

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

    params.solver_params.ls_maxstep = dd_nonlin_ls_maxstep;

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
        dd.set_to_remembered_solution();
        dd.set_simulation_voltage("gate", *it);
        sweep_drain(vds_stop, vds_steps, dd, *it);
      }

      vector<double>::iterator zero =
        find_if(voltages.begin(), first_positive,
            bind2nd(greater<double>(), -delta_v));

      if (zero != first_positive)
        sweep_drain(vds_stop, vds_steps, dd, 0.0);

      it = zero;
      if (it != voltages.begin())
      {
        bool restart = true;
        do
        {
          --it;
          if (!restart)
            dd.set_to_remembered_solution();
          dd.set_simulation_voltage("gate", *it);
          sweep_drain(vds_stop, vds_steps, dd, *it, restart);
          restart = false;
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
        const map<const ElectricalContact*, double>& curr =
          dd.get_boundary_currents();
        file << *it << " "
          << (*curr.find(dd.get_device().get_boundary("gate"))).second
          << " "
          << (*curr.find(dd.get_device().get_boundary("drain"))).second
          << " "
          << (*curr.find(dd.get_device().get_boundary("source"))).second
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
          const map<const ElectricalContact*, double>& curr =
            dd.get_boundary_currents();
          file << *it << " "
            << (*curr.find(dd.get_device().get_boundary("gate"))).second
            << " "
            << (*curr.find(dd.get_device().get_boundary("drain"))).second
            << " "
            << (*curr.find(dd.get_device().get_boundary("source"))).second
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
    DriftDiffusion& dd, double vg, bool restart)
{

  ostringstream filename;
  filename.precision(3);
  filename << "output/ids_" << fixed << vg << "V.dat";
  ofstream file;
  file.open(filename.str().c_str());

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

  vector<double> densities;
  vector<string> names;

  bool remember = true;
  vector<double>::iterator it = voltages.begin();
  for ( ; it != voltages.end(); ++it)
  {
    dd.set_simulation_voltage("drain", *it);
    cout << "Vgs = " << vg << " Vds = " << *it << "\n" << flush;
    dd.solve(restart);
    if (remember)
      dd.remember_current_solution();
    remember = false;
    restart = false;
    const map<const ElectricalContact*, double>& curr =
      dd.get_boundary_currents();
    file << *it << " "
         << (*curr.find(dd.get_device().get_boundary("gate"))).second
         << " "
         << (*curr.find(dd.get_device().get_boundary("drain"))).second
         << " "
         << (*curr.find(dd.get_device().get_boundary("source"))).second
         << "\n" << flush;

    ostringstream f;
    f.precision(3);
    f << "_" << fixed << vg << "V_" << fixed << *it << "V.gmv";
    dd.build_densities(densities, names);
    GMVIO(dd.get_mesh()).write_nodal_data("output/densities"+f.str(),
        densities, names);
    dd.build_band_edges(densities, names);
    GMVIO(dd.get_mesh()).write_nodal_data("output/band_edges"+f.str(),
        densities, names);
    dd.build_electric_field(densities, names);
    GMVIO_cell(dd.get_mesh()).write_ascii_cell_data("output/field"+f.str(),
        densities, names);
    dd.build_current_density(densities, names);
    GMVIO_cell(dd.get_mesh()).write_ascii_cell_data("output/current"+f.str(),
        densities, names);
  }
  
  file.close();
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
        bool found = true;
        AutoPtr<Elem> side = elem->build_side(s);
        for (int i = 0; i < side->n_nodes(); i++)
        {
          if (find(n_begin, n_end, side->node(i)) == n_end)
            found = false;
        }
        if (found)
          data.set_data(BoundaryData::ElementSide(elem, s), desc);
      }
    }
  }
}


