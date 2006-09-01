// the following _HAS_ to be included first
#include "Read_MSH.h"

#include "Dopant.h"
#include "RecombinationModelInterface.h"
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

    int mat_sys = input_file("material_system",0);
    double in_content = input_file("in_content", 0.15);
    double al_content = input_file("al_content", 0.48);
    double ga_content = input_file("ga_content", 0.47);
    int growth_dir = input_file("growth_dir",1);

    // drift-diffusion options
    GetPot dd_opt("dd.in");
    int characteristic = dd_opt("characteristic", 0);
    double n_doping = dd_opt("n_doping", 0.0);
    double p_doping = dd_opt("p_doping", 0.0);

    double start_voltage = dd_opt("vg_start", 0.0);
    double stop_voltage = dd_opt("vg_stop", 0.0);
    unsigned int voltage_steps = dd_opt("vg_steps", 1);

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

 
    mesh.print_info();


    /*****************************************************
     *
     * Setup of strain parameters
     *
     *****************************************************/

    bool periodicity[3];

    rotated_crystal cryst;

    if ((mat_sys != 0) && (mat_sys != 4))
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
      cryst.set_xyz_mil_direction("y", 1,  0, -1, 0) ;
      cryst.set_xyz_mil_direction("z", 1,  -2, 1, 0) ;
      cryst.set_xyz_mil_direction("x",  0,  0, 0,  1) ;
    }

    
    Macrostrain::strain_param a_strain;
    Macrostrain::strain_param b_strain;
    Macrostrain::strain_param c_strain;
    Piezoelectricity a_piezo;
    Piezoelectricity b_piezo;
    Piezoelectricity c_piezo;

    std::map<unsigned int, Macrostrain::strain_param*> strain_params;
    strain_params[1] = &a_strain;
    strain_params[2] = &b_strain;
    strain_params[3] = &c_strain;

    std::map<unsigned int, Piezoelectricity*> piezodata;
    piezodata[1] = &a_piezo;
    piezodata[2] = &b_piezo;
    piezodata[3] = &c_piezo;


    double x_al = al_content;
    double x_ga = ga_content;
    double x_in = 1 - al_content;
    
    if (mat_sys == 2)
    {
      // GaInAs - AlInAs
      cout << "Using GaInAs - AlInAs\n" << flush;

      a_strain.C_tensor.set_moduli(
          alloy(83.290, 125.0, x_in),
          alloy(45.260, 53.4, x_in),
          alloy(39.590, 54.2, x_in)); //InAlAs
      b_strain.C_tensor.set_moduli(
          alloy(122.1, 83.290, x_ga),
          alloy(56.6, 45.260, x_ga),
          alloy(60.0, 39.590, x_ga)); //GaInAs
      c_strain.C_tensor.set_moduli(
          alloy(83.290, 125.0, x_in),
          alloy(45.260, 53.4, x_in),
          alloy(39.590, 54.2, x_in)); //InAlAs

      a_piezo.set_moduli(alloy(-0.044, -0.015, x_in)); //InAs
      b_piezo.set_moduli(alloy(-0.16, -0.044, x_ga)); //GaAs  C/m^2
      c_piezo.set_moduli(alloy(-0.044, -0.015, x_in)); //InAs
      a_piezo.set_pyro_module(0.0);
      b_piezo.set_pyro_module(0.0);
      c_piezo.set_pyro_module(0.0);

      a_strain.crystal = cryst;
      b_strain.crystal = cryst;
      c_strain.crystal = cryst;

      std::vector<int> x_dir(3);
      std::vector<int> y_dir(3);
      // growth direction
      //x_dir[0] = growth_dir;  x_dir[1] = 1;  x_dir[2] =   1;
      //x_dir[0] = 1;  x_dir[1] = growth_dir;  x_dir[2] =   growth_dir;
      //y_dir[0] = 0;  y_dir[1] = 1;  y_dir[2] =  -1;
      y_dir[0] = 0;  y_dir[1] = 1;  y_dir[2] =  0;
      x_dir[0] = 1;  x_dir[1] = 0;  x_dir[2] =  0;

      a_strain.crystal.set_lat_const(alloy(6.05830, 5.6611, x_in));
      a_strain.crystal.calculate_lat_consts();
      a_strain.crystal.calculate_rot_matrix(x_dir, y_dir);
      a_strain.C_tensor.rotate_to_calc_system(a_strain.crystal.RotMatrix);
      
      b_strain.crystal.set_lat_const(alloy(5.65325, 6.05830, x_ga));
      b_strain.crystal.calculate_lat_consts();
      b_strain.crystal.calculate_rot_matrix(x_dir, y_dir);
      b_strain.C_tensor.rotate_to_calc_system(b_strain.crystal.RotMatrix);

      c_strain.crystal.set_lat_const(alloy(6.05830, 5.6611, x_in));
      c_strain.crystal.calculate_lat_consts();
      c_strain.crystal.calculate_rot_matrix(x_dir, y_dir);
      c_strain.C_tensor.rotate_to_calc_system(c_strain.crystal.RotMatrix);

    }
    else if (mat_sys == 1)
    {
      // GaAs - InAs
      cout << "Using GaAs - InAs\n" << flush;

      a_strain.C_tensor.set_moduli(122.1,  56.6,   60.0); //GaAs
      b_strain.C_tensor.set_moduli(83.290,  45.260,  39.590); //InAs
      c_strain.C_tensor.set_moduli(122.1,  56.6,   60.0); //GaAs

      a_piezo.set_moduli(-0.16); //GaAs  C/m^2
      b_piezo.set_moduli(-0.044); //InAs
      c_piezo.set_moduli(-0.16); //GaAs
      a_piezo.set_pyro_module(0.0);
      b_piezo.set_pyro_module(0.0);
      c_piezo.set_pyro_module(0.0);


      a_strain.crystal = cryst;
      b_strain.crystal = cryst;
      c_strain.crystal = cryst;

      std::vector<int> x_dir(3);
      std::vector<int> y_dir(3);
      // growth direction
      x_dir[0] = growth_dir;  x_dir[1] = 1;  x_dir[2] =   1;
      //x_dir[0] = 1;  x_dir[1] = growth_dir;  x_dir[2] =   growth_dir;
      y_dir[0] = 0;  y_dir[1] = 1;  y_dir[2] =  -1;

      a_strain.crystal.set_lat_const(5.65325);
      a_strain.crystal.calculate_lat_consts();
      a_strain.crystal.calculate_rot_matrix(x_dir, y_dir);
      a_strain.C_tensor.rotate_to_calc_system(a_strain.crystal.RotMatrix);

      b_strain.crystal.set_lat_const(6.05830);
      b_strain.crystal.calculate_lat_consts();
      b_strain.crystal.calculate_rot_matrix(x_dir, y_dir);
      b_strain.C_tensor.rotate_to_calc_system(b_strain.crystal.RotMatrix);

      c_strain.crystal.set_lat_const(5.65325);
      c_strain.crystal.calculate_lat_consts();
      c_strain.crystal.calculate_rot_matrix(x_dir, y_dir);
      c_strain.C_tensor.rotate_to_calc_system(c_strain.crystal.RotMatrix);
    }
    else if (mat_sys == 3)
    {
      // GaAs - InAs
      cout << "Using GaAs - InAs\n" << flush;

      a_strain.C_tensor.set_moduli(83.290,  45.260,  39.590); //InAs
      b_strain.C_tensor.set_moduli(122.1,  56.6,   60.0); //GaAs
      c_strain.C_tensor.set_moduli(83.290,  45.260,  39.590); //InAs

      a_piezo.set_moduli(-0.044); //InAs
      b_piezo.set_moduli(-0.16); //GaAs  C/m^2
      c_piezo.set_moduli(-0.044); //InAs
      a_piezo.set_pyro_module(0.0);
      b_piezo.set_pyro_module(0.0);
      c_piezo.set_pyro_module(0.0);


      a_strain.crystal = cryst;
      b_strain.crystal = cryst;
      c_strain.crystal = cryst;

      std::vector<int> x_dir(3);
      std::vector<int> y_dir(3);
      // growth direction
      x_dir[0] = 1;  x_dir[1] = 0;  x_dir[2] =   0;
      y_dir[0] = 0;  y_dir[1] = 1;  y_dir[2] =  0;

      a_strain.crystal.set_lat_const(6.05830);
      a_strain.crystal.calculate_lat_consts();
      a_strain.crystal.calculate_rot_matrix(x_dir, y_dir);
      a_strain.C_tensor.rotate_to_calc_system(a_strain.crystal.RotMatrix);

      b_strain.crystal.set_lat_const(5.65325);
      b_strain.crystal.calculate_lat_consts();
      b_strain.crystal.calculate_rot_matrix(x_dir, y_dir);
      b_strain.C_tensor.rotate_to_calc_system(b_strain.crystal.RotMatrix);

      c_strain.crystal.set_lat_const(6.05830);
      c_strain.crystal.calculate_lat_consts();
      c_strain.crystal.calculate_rot_matrix(x_dir, y_dir);
      c_strain.C_tensor.rotate_to_calc_system(c_strain.crystal.RotMatrix);
    }
    else
    {

      a_strain.crystal = cryst;
      b_strain.crystal = cryst;
      c_strain.crystal = cryst;

      std::vector<int> x_dir(4);
      std::vector<int> y_dir(4);

      y_dir[0] =  1;  y_dir[1] = 0;  y_dir[2] = -1;  y_dir[3] =  0;
      x_dir[0] =  0;  x_dir[1] = 0;  x_dir[2] =  0;  x_dir[3] =  1;


      // GaN
      a_strain.C_tensor.set_moduli(390.0, 145.0, 106.0, 398.0, 105.0);
      c_strain.C_tensor.set_moduli(390.0, 145.0, 106.0, 398.0, 105.0);
      //
      a_piezo.set_moduli(1.27, -0.35, -0.3); // C/m^2
      c_piezo.set_moduli(1.27, -0.35, -0.3); // C/m^2
      //
      a_piezo.set_pyro_module(-0.034);
      c_piezo.set_pyro_module(-0.034);
      //
      a_strain.crystal.set_lat_const(0.3189, 0.5185);
      c_strain.crystal.set_lat_const(0.3189, 0.5185);
      a_strain.crystal.calculate_lat_consts();
      c_strain.crystal.calculate_lat_consts();
      a_strain.crystal.calculate_rot_matrix(x_dir, y_dir); 
      c_strain.crystal.calculate_rot_matrix(x_dir, y_dir); 
      a_strain.C_tensor.rotate_to_calc_system(a_strain.crystal.RotMatrix);
      c_strain.C_tensor.rotate_to_calc_system(c_strain.crystal.RotMatrix);

      if (mat_sys == 4)
      {
        // AlN
        b_strain.C_tensor.set_moduli(396.0, 137.0, 108.0, 373.0, 116.0);
        b_piezo.set_moduli(1.79, -0.5, -0.48);
        b_piezo.set_pyro_module(-0.09);
        b_strain.crystal.set_lat_const(0.3112, 0.4982);
        b_strain.crystal.calculate_lat_consts();
        b_strain.crystal.calculate_rot_matrix(x_dir, y_dir);
        b_strain.C_tensor.rotate_to_calc_system(b_strain.crystal.RotMatrix);
      }
      else
      {    
        // InGaN
        // first value: InN, second value: GaN
        b_strain.C_tensor.set_moduli(
            alloy(223.0, 390.0, in_content),
            alloy(115.0, 145.0, in_content),
            alloy(92.0, 106.0, in_content),
            alloy(224.0, 398.0, in_content),
            alloy(48.0, 105.0, in_content));

        b_piezo.set_moduli(
            alloy(0.97, 1.27, in_content),
            alloy(-0.57, -0.35, in_content),
            alloy(-0.4, -0.3, in_content));
        b_piezo.set_pyro_module(
            alloy(-0.042, -0.034, in_content, -0.037));

        b_strain.crystal.set_lat_const(
            alloy(0.3545, 0.3189, in_content),
            alloy(0.5703, 0.5185, in_content));
        b_strain.crystal.calculate_lat_consts();
        b_strain.crystal.calculate_rot_matrix(x_dir, y_dir);
        b_strain.C_tensor.rotate_to_calc_system(b_strain.crystal.RotMatrix);
      }
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

    double stress_value = strain_opt( "stress",0.0);
    map<unsigned int, double> stress_map;
    stress_map[1] = stress_value;


    EquationSystems eqsys(mesh);
    Macrostrain strain_calculation(opt, eqsys, "strainsys");

    strain_calculation.define_substrate_bc(2);
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

    StrainedSemiconductorModel n_gan(&strain_calculation);    
    StrainedSemiconductorModel p_gan(&strain_calculation);
    if (mat_sys == 1)
    {
      n_gan.set_data_file("materials/GaAs.dat");
      p_gan.set_data_file("materials/GaAs.dat");
    }
    else if (mat_sys == 2)
    {
      n_gan.set_data_file("materials/InAs.dat");
      p_gan.set_data_file("materials/InAs.dat");
    }
    else if (mat_sys == 3)
    {
      n_gan.set_data_file("materials/InAs.dat");
      p_gan.set_data_file("materials/InAs.dat");
    }
    else
    {
      n_gan.set_data_file("materials/GaN.dat");
      p_gan.set_data_file("materials/GaN.dat");
    }

    
    ModelOptions opts;
    opts["tau_n"] = "1e-10";
    opts["tau_p"] = "1e-10";
    opts["C"] = "1e-10";
    RecombinationModelInterface* rm =
      RecombinationModelInterface::create("SRH", opts);
    n_gan.add_recombination_model(rm);
    rm = RecombinationModelInterface::create("direct", opts);
    n_gan.add_recombination_model(rm);

    rm = RecombinationModelInterface::create("SRH", opts);
    p_gan.add_recombination_model(rm);
    rm = RecombinationModelInterface::create("direct", opts);
    p_gan.add_recombination_model(rm);
    
    // InGaN
    StrainedSemiconductorModel ingan(&strain_calculation);
    if (mat_sys == 1)
      ingan.set_data_file("materials/InAs.dat");
    else if (mat_sys == 3)
      ingan.set_data_file("materials/GaAs.dat");
    else if (mat_sys == 4)
      ingan.set_data_file("materials/AlN.dat");

    opts["tau_n"] = "1e-7";
    opts["tau_p"] = "1e-7";
    opts["C"] = "1e-10";
    rm = RecombinationModelInterface::create("SRH", opts);
    ingan.add_recombination_model(rm);
    rm = RecombinationModelInterface::create("direct", opts);
    ingan.add_recombination_model(rm);

    n_gan.add_dopant(new Dopant(n_doping, 0.025, 2, Dopant::N_TYPE));
    p_gan.add_dopant(new Dopant(p_doping, 0.01, 4, Dopant::P_TYPE));

    n_gan.set_statistics(TiberCad::FERMIDIRAC);
    p_gan.set_statistics(TiberCad::FERMIDIRAC);
    ingan.set_statistics(TiberCad::FERMIDIRAC);

    Dummy d;
    n_gan.read_database(d);
    p_gan.read_database(d);
    ingan.read_database(d);

    if (mat_sys == 2)
    {
      n_gan.build_alloy("materials/AlAs.dat",
          "materials/AlInAs_bow.dat", x_al);
      p_gan.build_alloy("materials/AlAs.dat",
          "materials/AlInAs_bow.dat", x_al);
      ingan.build_alloy("materials/GaAs.dat",
          "materials/GaInAs_bow.dat", x_ga);
    }
    else if (mat_sys == 0)
      ingan.build_alloy("materials/InN.dat",
          "materials/InGaN_bow.dat", in_content);
    
    n_gan.set_mobilities(200, 50);
    p_gan.set_mobilities(200, 50);
    ingan.set_mobilities(200, 50);

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



    cout << "Solving strain... \n" << flush;
    strain_calculation.solve();
    strain_calculation.output_piezo("Piezo.gmv");


    OhmicContact anode("anode");
    anode.set_zero_derivative_bc(POTENTIAL);
    anode.set_zero_derivative_bc(FERMIE);
    OhmicContact cathode("cathode");
    cathode.set_zero_derivative_bc(POTENTIAL);
    cathode.set_zero_derivative_bc(FERMIH);


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
            set_boundary(boundary_data, nodes, &anode, mesh);
            break;
          case 2:
            set_boundary(boundary_data, nodes, &cathode, mesh);
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
    dd.set_equation_systems(&eqsys);
    dd.init();

    DriftDiffusion::Options& params = dd.get_options();
    params.max_refinement_steps = dd_max_r_steps;
    params.solver_params.nonlinear_max_iterations = 100;
    params.solver_params.linear_max_iterations = dd_lin_max_it;
    params.solver_params.nonlinear_tolerance = dd_nonlin_rtol;
    params.solver_params.nonlinear_abs_tolerance = dd_nonlin_atol;
    params.solver_params.ls_maxstep = 0.05;
    params.solver_params.linear_tolerance = dd_lin_rtol;
    params.solver_params.linear_abs_tolerance = dd_lin_atol;
    params.integration_order =
      static_cast<libMeshEnums::Order>(integration_order);

    params.solver_params.ksp_type = KSPGMRES;
    params.solver_params.pc_type = PCILU;
    //params.artificial_drift = true;
    params.local_scaling = true;
    

    // mesh drawn in um
    params.mesh_units = mesh_units;

    dd.enable_mesh_refinement();

    /*****************************************************/

    


    cout << "Solving drift-diffusion... \n" << flush;

    dd.set_simulation_voltage("cathode", 0.0);
    dd.set_simulation_voltage("anode", 0.0);
    
    cout << "Solving equilibrium... " << flush;
    params.coupling = POISSON;
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
    cout << "Barrier:\n";
    n_gan.print_info();
    cout << "\nWell:\n";
    ingan.print_info();

    params.solver_params.nonlinear_max_iterations = dd_nonlin_max_it;
    params.solver_params.ls_maxstep = dd_nonlin_ls_maxstep;

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

