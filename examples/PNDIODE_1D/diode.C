#include "ElementData.h"
#include "OhmicContact.h"
#include "BoundaryData.h"
#include "DDevice.h"
#include "DriftDiffusion.h"
#include "SemiconductorModel.h"

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

void set_boundary(BoundaryData& data, double x,
    ElectricalContact* desc, const Mesh& mesh);


class Dummy {};


int main (int argc, char** argv)
{
  libMesh::init(argc, argv);
  {

    double fac = 1e-0;
    
    GetPot input_file("diode.in");

    const string approx_order = input_file("approximation_order", "FIRST");

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

    vector<double> phys_reg_ID(1);
    phys_reg_ID[0] = 101; // p
    //phys_reg_ID[1] = 102; // n

    
    Mesh mesh(1);
    //MeshData_elements meshdata(mesh);
    //meshdata.enable_compatibility_mode();
    //readmesh.read_mesh_and_data(mesh, meshdata);

    if (approx_order == "FIRST")
    {
      MeshTools::Generation::build_line(mesh, 
          initial_density,
          xmin, xmax,
          EDGE2);
    }
    else
    {
      MeshTools::Generation::build_line(mesh, 
          initial_density,
          xmin, xmax,
          EDGE3);
    }
    mesh.print_info();


    SemiconductorModel nside;    
    nside.set_data_file("Si.dat");
    if (statistics == "FD")
      nside.set_statistics(TiberCad::FERMIDIRAC);
    else
      nside.set_statistics(TiberCad::BOLTZMANN);

    nside.set_mobilities(800, 200);
    nside.add_recombination_model(SRH);
    nside.set_SRH_parameters(1e-9, 1e-9);

    SemiconductorModel pside(nside);

    nside.set_n_dopant(Dopant(n_doping, 0.025, 2));
    pside.set_p_dopant(Dopant(p_doping, 0.01, 4));

    Dummy d;
    nside.read_database(d);
    pside.read_database(d);


    ElementData element_data;

    {
      Mesh::const_element_iterator it = mesh.elements_begin();
      const Mesh::const_element_iterator end = mesh.elements_end();

      for ( ; it != end; ++it)
      {
        const Elem* elem = *it;

        Point c = elem->centroid();
        if (c(0) < 0.0)
          element_data.set_data(elem, &pside);
        else
          element_data.set_data(elem, &nside);
      }

    }

    OhmicContact anode("anode");
    anode.set_zero_derivative_bc(FERMIE);
    OhmicContact cathode("cathode");
    cathode.set_zero_derivative_bc(FERMIH);

    BoundaryData boundary_data;
    set_boundary(boundary_data, xmin, &anode, mesh);
    set_boundary(boundary_data, xmax, &cathode, mesh);



    DD::Device device(&mesh, &element_data, &boundary_data);
    bool device_integrity = device.check_integrity();
    if (device_integrity)
      cout << "Device ok.\n\n";
    else
    {
      cout << "Device bad.\n\n";
      return 1;
    }

    MeshRefinement ref(mesh);
    /*
    for (int i = 0; i < 1; i++)
    {
      {
        MeshBase::element_iterator it = mesh.elements_begin();
        const MeshBase::element_iterator end = mesh.elements_end();
        for ( ; it != end; ++it)
        {
          Elem* elem = *it;

          double y = elem->centroid()(0);
          if ((y > -2.0) && (y < 2.0))
            elem->set_refinement_flag(Elem::REFINE);
        }
        ref.refine_elements();
      }
    }
    for (int i = 0; i < 1; i++)
    {
      {
        MeshBase::element_iterator it = mesh.elements_begin();
        const MeshBase::element_iterator end = mesh.elements_end();
        for ( ; it != end; ++it)
        {
          Elem* elem = *it;

          double y = elem->centroid()(0);
          if ((y > -1.0) && (y < 1.0))
            elem->set_refinement_flag(Elem::REFINE);
        }
        ref.refine_elements();
      }
    }
    */
    for (int i = 0; i < 1; i++)
    {
      {
        MeshBase::element_iterator it = mesh.elements_begin();
        const MeshBase::element_iterator end = mesh.elements_end();
        for ( ; it != end; ++it)
        {
          Elem* elem = *it;

          double y = elem->centroid()(0);
          if ((y > -0.5) && (y < 0.5))
            elem->set_refinement_flag(Elem::REFINE);
        }
        ref.refine_elements();
      }
    }
    for (int i = 0; i < 1; i++)
    {
      {
        MeshBase::element_iterator it = mesh.elements_begin();
        const MeshBase::element_iterator end = mesh.elements_end();
        for ( ; it != end; ++it)
        {
          Elem* elem = *it;

          double y = elem->centroid()(0);
          if ((y > -0.2) && (y < 0.2))
            elem->set_refinement_flag(Elem::REFINE);
        }
        ref.refine_elements();
      }
    }
    /*
    for (int i = 0; i < 1; i++)
    {
      {
        MeshBase::element_iterator it = mesh.elements_begin();
        const MeshBase::element_iterator end = mesh.elements_end();
        for ( ; it != end; ++it)
        {
          Elem* elem = *it;

          double y = elem->centroid()(0);
          if ((y > -0.05) && (y < 0.05))
            elem->set_refinement_flag(Elem::REFINE);
        }
        ref.refine_elements();
      }
    }
    for (int i = 0; i < 1; i++)
    {
      {
        MeshBase::element_iterator it = mesh.elements_begin();
        const MeshBase::element_iterator end = mesh.elements_end();
        for ( ; it != end; ++it)
        {
          Elem* elem = *it;

          double y = elem->centroid()(0);
          if ((y > -0.005) && (y < 0.005))
            elem->set_refinement_flag(Elem::REFINE);
        }
        ref.refine_elements();
      }
    }
    for (int i = 0; i < 1; i++)
    {
      {
        MeshBase::element_iterator it = mesh.elements_begin();
        const MeshBase::element_iterator end = mesh.elements_end();
        for ( ; it != end; ++it)
        {
          Elem* elem = *it;

          double y = elem->centroid()(0);
          if ((y > -0.002) && (y < 0.002))
            elem->set_refinement_flag(Elem::REFINE);
        }
        ref.refine_elements();
      }
    }
    */


  
    EquationSystems eqsys(mesh);
    DriftDiffusion dd(&device);
    dd.set_equation_systems(&eqsys);

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


    if (approx_order == "FIRST")
      params.approximation_order = FIRST;
    else if (approx_order == "SECOND")
      params.approximation_order = SECOND;

    // mesh drawn in um
    params.mesh_units = mesh_units;

    dd.enable_mesh_refinement();

    
    dd.set_simulation_voltage("cathode", 0.0);
    dd.set_simulation_voltage("anode", 0.0);
    {
      vector<double> densities;
      vector<string> names;
      cout << "Solving equilibrium... " << flush;
      dd.solve();
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

    const Scaling& sc = dd.get_scaling();
    cout << "Scaling parameters:\n";
    cout << "     phi0: " << sc.get_potential_scaling() << "\n";
    cout << "     x0  : " << sc.get_length_scaling() << "\n";
    cout << "     mu0 : " << sc.get_mobility_scaling() << "\n";
    cout << "     C0  : " << sc.get_density_scaling() << "\n";
    cout << "     t0  : " << sc.get_time_scaling() << "\n";
    cout << "     R0  : " << sc.get_recombination_scaling() << "\n\n";



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
    vector<double> voltages(voltage_steps + 1);
    for (int i = 0; i <= voltage_steps; i++)
    {
      voltages[i] = start_voltage + i * step;
    }

    vector<double>::iterator first_positive =
      find_if(voltages.begin(), voltages.end(),
          bind2nd(greater<double>(), delta_v));

    map<double, vector<double> > iv_char;

    vector<double>::iterator it = first_positive;
    for ( ; it != voltages.end(); ++it)
    {
      dd.set_simulation_voltage("anode", *it);
      cout << " Solving U = " << *it << " V ... " << flush;
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
      iv_char[*it] = vector<double>(3);
      iv_char[*it][0] = (*curr.find(&cathode)).second;
      iv_char[*it][1] = (*curr.find(&anode)).second;
      iv_char[*it][2] = dd.get_artificial_boundary_current();
      cerr << "    I = " << (iv_char[*it][1] * fac) << " A/cm^2\n";
    }

    vector<double>::iterator zero =
      find_if(voltages.begin(), first_positive,
          bind2nd(greater<double>(), -delta_v));

    if (zero != first_positive)
      iv_char[0.0] = vector<double>(3, 0.0);

    it = zero;
    if (it != voltages.begin())
    {
      bool restart = true;
      do
      {
        --it;
        dd.set_simulation_voltage("anode", *it);
        cout << " Solving U = " << *it << " V ... " << flush;
        dd.enable_mesh_refinement();
        dd.solve(restart);
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
        iv_char[*it] = vector<double>(3);
        iv_char[*it][0] = (*curr.find(&cathode)).second;
        iv_char[*it][1] = (*curr.find(&anode)).second;
        iv_char[*it][2] = dd.get_artificial_boundary_current();
        cerr << "    I = " << (iv_char[*it][1] * fac) << " A\n";
      }
      while (it != voltages.begin());
    }
    
    ofstream file;
    file.open("output/iv_char.dat");
    //file << "# V      A/cm\n";
    map<double, vector<double> >::const_iterator iv_it = iv_char.begin();
    for ( ; iv_it != iv_char.end(); ++iv_it)
    {
      //cout << (*iv_it).first << "V - "
      //     << (*iv_it).second << " A/cm\n";
      const vector<double>& curr = (*iv_it).second;
      file << (*iv_it).first << "  "
           << curr[0] << "  "
           << curr[1] << "  "
           << curr[2] << "\n";
    }

  }

  return libMesh::close();
}



void set_boundary(BoundaryData& data, double x,
    ElectricalContact* desc, const Mesh& mesh)
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
          data.set_data(BoundaryData::ElementSide(elem, n), desc);
        }
      }
    }

    ++el;
  }
}

