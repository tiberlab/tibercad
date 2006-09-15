// the following _HAS_ to be included first
#include "Read_MSH.h"

#include "ElementData.h"
#include "OhmicContact.h"
#include "BoundaryData.h"
#include "DDevice.h"
#include "DriftDiffusion.h"
#include "Dopant.h"
#include "StrainedSemiconductorModel.h"
#include "ExcitonTransport.h"
#include "ExcitonModel.h"
#include "RecombinationModelInterface.h"
#include "ExcitonDissociation.h"
#include "macrostrain.h"
#include "WzDDsemiconductor.h"

#include "mesh.h"
#include "mesh_modification.h"
#include "mesh_data.h"
#include "mesh_generation.h"
#include "elem.h"
#include "getpot.h"
#include "gnuplot_io.h"
#include "equation_systems.h"

#include <algorithm>

using namespace std;
using namespace DriftDiffusionDefs;

void set_boundary(BoundaryData& data, const vector<unsigned int>& nodes,
    ElectricalContact* desc, const Mesh& mesh);

class Dummy {};


double alloy(double a, double b, double xa, double bowing = 0.0)
{
  return b + (a - b) * xa - bowing * xa *  (1 - xa);
};

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
    double mu_e = input_file("electron_mobility", 800);
    double mu_h = input_file("hole_mobility", 200);
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

    double in_content = input_file("in_content", 0.14);


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



    vector<unsigned int> phys_reg_ID(3);
    phys_reg_ID[0] = 1; // n
    phys_reg_ID[1] = 2; // QW
    phys_reg_ID[2] = 3; // p
    
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
    
    mesh.print_info();

    EquationSystems eqsys(mesh);

    /*****************************************************
     *
     * Setup of strain parameters
     *
     *****************************************************/

    bool periodicity[3];

    rotated_crystal cryst;
    // AlGaN - GaN
    cryst.set_cryst_type("hex");
    cryst.set_xyz_mil_direction("y", 1,  0, -1, 0) ;
    cryst.set_xyz_mil_direction("z", 1,  -2, 1, 0) ;
    cryst.set_xyz_mil_direction("x",  0,  0, 0,  1) ;
    
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


    Macrostrain strain_calculation(opt, eqsys, "strainsys");

    strain_calculation.define_substrate_bc(1);
    strain_calculation.define_BC_map(boundary_nodes);
    strain_calculation.define_stress_value(stress_map);

    strain_calculation.assign_mesh_data(meshdata);
    strain_calculation.define_strain_parameters(strain_params);
    strain_calculation.define_piezo_moduli(piezodata);



    // Solve strain
    strain_calculation.solve();
    strain_calculation.output_piezo("output/piezo.gmv");
    strain_calculation.output_strain("output/strain.gmv");



    /****************************************************************
     *
     * Drift diffusion definitions
     *
     ****************************************************************/


    Dummy d;

    StrainedSemiconductorModel nside(&strain_calculation);    
    nside.set_data_file(material);
    nside.read_database(d);

    nside.add_dopant(new Dopant(n_doping, 0.025, 2, Dopant::N_TYPE));
    nside.set_mobilities(mu_e, mu_h);

    ModelOptions opts;
    opts["tau_n"] = tau_n;
    opts["tau_p"] = tau_n;
    opts["C"] = C_direct;
    RecombinationModelInterface* rm =
      RecombinationModelInterface::create("SRH", opts);
    nside.add_recombination_model(rm);
    rm = RecombinationModelInterface::create("direct", opts);
    nside.add_recombination_model(rm);
    opts["C"] = X_gen;
    rm = RecombinationModelInterface::create("exciton_generation", opts);
    nside.add_recombination_model(rm);


    StrainedSemiconductorModel pside(&strain_calculation);
    pside.set_data_file(material);
    pside.read_database(d);

    pside.add_dopant(new Dopant(5e16, 0.025, 2, Dopant::N_TYPE));
    pside.add_dopant(new Dopant(p_doping, 0.17, 4, Dopant::P_TYPE));
    pside.set_mobilities(mu_e, mu_h);

    opts["tau_n"] = tau_p;
    opts["tau_p"] = tau_p;
    opts["C"] = C_direct;
    rm = RecombinationModelInterface::create("SRH", opts);
    pside.add_recombination_model(rm);
    rm = RecombinationModelInterface::create("direct", opts);
    pside.add_recombination_model(rm);
    opts["C"] = X_gen;
    rm = RecombinationModelInterface::create("exciton_generation", opts);
    pside.add_recombination_model(rm);


    StrainedSemiconductorModel qw(&strain_calculation);
    qw.set_data_file(material);
    qw.read_database(d);

    qw.build_alloy("../material/InN.dat",
        "../material/InGaN_bow.dat", in_content);

    opts["tau_n"] = "1e-9";
    opts["tau_p"] = "1e-9";
    opts["C"] = C_direct;
    rm = RecombinationModelInterface::create("SRH", opts);
    qw.add_recombination_model(rm);
    rm = RecombinationModelInterface::create("direct", opts);
    qw.add_recombination_model(rm);
    opts["C"] = X_gen;
    rm = RecombinationModelInterface::create("exciton_generation", opts);
    qw.add_recombination_model(rm);

    qw.set_mobilities(mu_e, mu_h);


    if (statistics == "FD")
    {
      nside.set_statistics(TiberCad::FERMIDIRAC);
      pside.set_statistics(TiberCad::FERMIDIRAC);
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
          case 2:
            element_data.set_data(elem, &qw);
            break;
          case 3:
            element_data.set_data(elem, &pside);
            break;
          default:
            element_data.set_data(elem, &nside);
            break;
        }
      }
    }

    OhmicContact anode("anode");
    //SchottkyContact anode("anode");
    //anode.set_schottky_barrier_height(0.8);
    anode.set_zero_derivative_bc(FERMIE);
    OhmicContact cathode("cathode");
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

        if (it->first == 1)
          set_boundary(boundary_data, nodes, &anode, mesh);
        else
          set_boundary(boundary_data, nodes, &cathode, mesh);
      }
    }


    ExcitonModel xpart;
    xpart.set_exciton_generation_model("exciton_generation");
    // this is needed, but a momentary quirk
    xpart.read_database(d);
    xpart.set_recombination_times(rec_rad, rec_nonrad);
    xpart.set_binding_energy(20.4e-3);
    xpart.set_effective_mass(1.02);  // (*me) in exciton model
    xpart.set_mobility(1500);




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
    dd.set_equation_systems(&eqsys);
    dd.init();

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


    ExcitonTransport ex(&device);
    ex.set_equation_systems(&eqsys);
    ex.init();
    ex.set_exciton_model(&xpart);
    ex.set_driftdiffusion(&dd);
    ex.set_initial_guess(1.5);

    ExcitonTransport::Options& exparams = ex.get_options();
    exparams.max_refinement_steps = refinement_steps;
    exparams.solver_params.nonlinear_max_iterations = 50;
    exparams.solver_params.linear_max_iterations = lin_max_it;
    exparams.solver_params.nonlinear_tolerance = 1e-19;
    exparams.solver_params.nonlinear_abs_tolerance = 1e-15;
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
      dd.solve();

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
      dd.build_electric_field(densities, names);
      GMVIO_cell(dd.get_mesh()).write_ascii_cell_data("output/field_eq",
            densities, names);

      ex.build_densities(densities, names);
      GnuPlotIO(dd.get_mesh()).write_nodal_data("output/excitons_eq",
            densities, names);
    }
    cout << "GaN material info:\n";
    cout << "n-side:\n";
    nside.print_info();
    cout << "p-side:\n";
    pside.print_info();

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
    opts["trapping_probability"] = trapping;
    RecombinationModelInterface* rma =
      RecombinationModelInterface::create("exciton_dissociation", opts);
    (static_cast<ExcitonDissociation*>(rma))->set_exciton_transport(&ex);
    RecombinationModelInterface* rmb =
      RecombinationModelInterface::create("exciton_dissociation", opts);
    (static_cast<ExcitonDissociation*>(rmb))->set_exciton_transport(&ex);

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
    

    bool done = false;
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

      if (it == first_positive)
      {
        //nside.add_recombination_model(rma);
        //pside.add_recombination_model(rmb);
      }


      if (*it >= 0.5)
      {
        if (!done)
        {
          nside.add_recombination_model(rma);
          pside.add_recombination_model(rmb);
          done = true;
        }

        try { ex.solve(); }
        catch (...) {}

        dd.solve();
        try { ex.solve(); }
        catch (...) {}
        {
          damping = pow(damping, 0.5);
          ostringstream s;
          s << damping;
          opts["damping"] = s.str();
          rma->set_model_options(opts);
          rmb->set_model_options(opts);
        }
        dd.solve();
        try { ex.solve(); }
        catch (...) {}
        /*
        {
          damping = pow(damping, 0.5);
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
      

      const map<const ElectricalContact*, double>& curr =
        dd.get_boundary_currents();
      file << *it << "  "
           << (*curr.find(&cathode)).second << "  "
           << (*curr.find(&anode)).second << "  "
           << dd.get_artificial_boundary_current() << "\n" << flush;
      cerr << "    I = " << (*curr.find(&cathode)).second << " A/cm\n";
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

        const map<const ElectricalContact*, double>& curr =
          dd.get_boundary_currents();
        file << *it << "  "
          << (*curr.find(&cathode)).second << "  "
          << (*curr.find(&anode)).second << "  "
          << dd.get_artificial_boundary_current() << "\n" << flush;
        cerr << "    I = " << (*curr.find(&cathode)).second << " A/cm\n";
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

