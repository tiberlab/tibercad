// the following _HAS_ to be included first
#include "Read_MSH.h"

#include "ElementData.h"
#include "OhmicContact.h"
#include "BoundaryData.h"
#include "DDevice.h"
#include "Dopant.h"
#include "DriftDiffusion.h"
#include "SemiconductorModel.h"
#include "RecombinationModelInterface.h"

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

void set_boundary(BoundaryData& data, const vector<unsigned int>& nodes,
    ElectricalContact* desc, const Mesh& mesh);



class Dummy {};


int main (int argc, char** argv)
{
  libMesh::init(argc, argv);
  {

    GetPot input_file("diode.in");

    string meshfile = input_file("meshfile", "");

    double start_voltage = input_file("start_voltage", 0.0);
    double stop_voltage = input_file("stop_voltage", 0.0);
    unsigned int voltage_steps = input_file("voltage_steps", 1);


    const string method = input_file("simulation_method", "NEWTON");
    const string statistics = input_file("statistics", "B");

    double nonlin_rtol = input_file("nonlinear_tolerance", 1e-9);
    double lin_rtol = input_file("linear_tolerance", 1e-12);
    int integration_order = input_file("integration_order", 5);
    int nonlin_max_it = input_file("nonlinear_max_it", 15);
    int lin_max_it = input_file("linear_max_it", 500);
    double nonlin_ls_maxstep = input_file("nonlinear_ls_maxstep", 0.025);

    const double mesh_units = input_file("mesh_units", 1e-4);
    double xmin = input_file("xmin", -0.5);
    double xmax = input_file("xmax", 0.5);

    double n_doping = input_file("n_doping", 1e18);
    double p_doping = input_file("p_doping", 1e18);

    unsigned int initial_density = input_file("initial_density", 20);

    unsigned int refinement_steps = input_file("max_refinement_steps", 0);
    double refine_frac = input_file("refine_fraction", 0.7);
    double coarsen_frac = input_file("coarsen_fraction", 0.3);

    vector<unsigned int> phys_reg_ID(2);
    phys_reg_ID[0] = 1; // p
    phys_reg_ID[1] = 2; // n
    
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
    

    mesh.print_info();

    Dummy d;

    SemiconductorModel nside;    
    nside.set_data_file("Si.dat");
    nside.read_database(d);

    nside.add_dopant(new Dopant(n_doping, 0.025, 2, Dopant::N_TYPE));
    nside.set_mobilities(283.8, 160.3);

    ModelOptions opts;
    opts["tau_n"] = "1e-7";
    opts["tau_p"] = "3e-8";
    RecombinationModelInterface* rm =
      RecombinationModelInterface::create("SRH", opts);
    nside.add_recombination_model(rm);


    SemiconductorModel pside;
    pside.set_data_file("Si.dat");
    pside.read_database(d);

    pside.add_dopant(new Dopant(p_doping, 0.01, 4, Dopant::P_TYPE));
    pside.set_mobilities(283.8, 160.3);

    opts["tau_n"] = "1e-7";
    opts["tau_p"] = "3e-8";
    rm = RecombinationModelInterface::create("SRH", opts);
    pside.add_recombination_model(rm);

    if (statistics == "FD")
    {
      nside.set_statistics(TiberCad::FERMIDIRAC);
      pside.set_statistics(TiberCad::FERMIDIRAC);
    }
    else
    {
      nside.set_statistics(TiberCad::BOLTZMANN);
      pside.set_statistics(TiberCad::BOLTZMANN);
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
            element_data.set_data(elem, &pside);
            break;
          default:
            element_data.set_data(elem, &nside);
            break;
        }
      }
    }

    OhmicContact anode("anode");
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


    DD::Device device(&mesh, &element_data, &boundary_data);
    bool device_integrity = device.check_integrity();
    if (device_integrity)
      cout << "Device ok.\n\n";
    else
    {
      cout << "Device bad.\n\n";
      return 1;
    }

  
    EquationSystems eqsys(mesh);
    DriftDiffusion dd(&device);
    dd.set_equation_systems(&eqsys);
    dd.init();

    DriftDiffusion::Options& params = dd.get_options();
    params.max_refinement_steps = refinement_steps;
    params.solver_params.nonlinear_max_iterations = nonlin_max_it;
    params.solver_params.linear_max_iterations = lin_max_it;
    params.solver_params.nonlinear_tolerance = nonlin_rtol;
    params.solver_params.ls_maxstep = nonlin_ls_maxstep;
    params.solver_params.linear_tolerance = lin_rtol;
    params.integration_order =
      static_cast<libMeshEnums::Order>(integration_order);
    params.refine_fraction = refine_frac;
    params.coarsen_fraction = coarsen_frac;
    params.solver_params.pc_type = PCILU;

    if (method == "GUMMEL")
    {
      params.solver_method = DriftDiffusion::GUMMEL;
      params.max_gummel_iterations = 2;
    }
    else
      params.solver_method = DriftDiffusion::NEWTON;


    // mesh drawn in um
    params.mesh_units = mesh_units;

    dd.enable_mesh_refinement();

    
    dd.set_simulation_voltage("cathode", 0.0);
    dd.set_simulation_voltage("anode", 0.0);
    {
      vector<double> densities;
      vector<string> names;
      cout << "Solving equilibrium...\n" << flush;
      dd.guess_equilibrium();
      dd.solve();
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
    cout << "Si properties:\n";
    nside.print_info();

    const Scaling& sc = dd.get_scaling();
    cout << "Scaling parameters:\n";
    cout << "     phi0: " << sc.get_potential_scaling() << "\n";
    cout << "     x0  : " << sc.get_length_scaling() << "\n";
    cout << "     mu0 : " << sc.get_mobility_scaling() << "\n";
    cout << "     C0  : " << sc.get_density_scaling() << "\n";
    cout << "     t0  : " << sc.get_time_scaling() << "\n";
    cout << "     R0  : " << sc.get_recombination_scaling() << "\n\n";



    cout << "\nBegin sweep...\n" << flush;
    params.coupling = FULLYCOUPLED;
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



