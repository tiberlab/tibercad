// the following _HAS_ to be included first
#include "Read_MSH.h"
#include "ReadISEGrid.h"

#include "Dopant.h"
#include "ElementData.h"
#include "OhmicContact.h"
#include "SchottkyContact.h"
#include "BoundaryData.h"
#include "Dopant.h"
#include "DDevice.h"
#include "DriftDiffusion.h"
#include "SemiconductorModel.h"
#include "RecombinationModelInterface.h"

#include "mesh.h"
#include "mesh_modification.h"
#include "mesh_data.h"
#include "mesh_generation.h"
#include "elem.h"
#include "getpot.h"
#include "gmv_io.h"
#include "tecplot_io.h"
#include "GMVIO_cell.h"
#include "equation_systems.h"

#include <algorithm>

using namespace std;
using namespace DriftDiffusionDefs;

void set_boundary(BoundaryData& data, const vector<unsigned int>& nodes,
    ElectricalContact* desc, const Mesh& mesh);

void sweep_drain(double stop, int steps, DriftDiffusion& dd, double vg);

class Dummy {};


int main (int argc, char** argv)
{
  libMesh::init(argc, argv);
  {

    GetPot input("mesfet.in");

    string meshfile = input("meshfile", "");
    string format = input("format", "gmsh");
    const double mesh_units = input("mesh_units", 1e-4);

    int fully = input("fully_coupled", 1);

    double schottky_barrier = input("schottky_barrier", 0.8);
    
    double vds_stop = input("vds_stop", 0.0);
    unsigned int vds_steps = input("vds_steps", 1);
    double vg_start = input("vg_start", 0.0);
    double vg_stop = input("vg_stop", 0.0);
    unsigned int vg_steps = input("vg_steps", 1);

    double vds = input("vds", 0.0);

    int characteristic = input("characteristic", 0);


    const string method = input("simulation_method", "NEWTON");
    const string statistics = input("statistics", "B");

    double nonlin_rtol = input("nonlinear_tolerance", 1e-9);
    double lin_rtol = input("linear_tolerance", 1e-12);
    int integration_order = input("integration_order", 5);
    int nonlin_max_it = input("nonlinear_max_it", 15);
    int lin_max_it = input("linear_max_it", 500);
    double nonlin_ls_maxstep =
      input("nonlinear_ls_maxstep", 0.025);


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
      ReadISEGrid readmesh(meshfile.c_str(), mesh, meshdata);
      readmesh.get_BC_data(boundary_nodes);
    }
    else
    {
      Read_MSH readmesh(meshfile, phys_reg_ID, BC_reg_ID, dim, mesh, meshdata);
      readmesh.get_BC_data(boundary_nodes);
    }

    mesh.print_info();

    Dummy d;

    SemiconductorModel sub;    
    sub.set_data_file("Si.dat");
    sub.read_database(d);
    {
      sub.add_dopant(new Dopant(1e17, 0.01, 4, Dopant::P_TYPE));
      //sub.set_mobilities(1416, 470);
      sub.set_mobilities(560, 280);
      ModelOptions opts;
      opts["tau_n"] = "1e-6";
      opts["tau_p"] = "3e-7";
      RecombinationModelInterface* rm =
        RecombinationModelInterface::create("SRH", opts);
      sub.add_recombination_model(rm);
    }



    SemiconductorModel schottky;
    schottky.set_data_file("Si.dat");
    schottky.read_database(d);
    { 
      schottky.add_dopant(new Dopant(1e18, 0.025, 2, Dopant::N_TYPE));
      schottky.set_mobilities(283.8, 160.3);
      ModelOptions opts;
      opts["tau_n"] = "1e-7";
      opts["tau_p"] = "3e-8";
      RecombinationModelInterface* rm =
        RecombinationModelInterface::create("SRH", opts);
      schottky.add_recombination_model(rm);
    }


    SemiconductorModel contact;
    contact.set_data_file("Si.dat");
    contact.read_database(d);
    {
      contact.add_dopant(new Dopant(5e19, 0.025, 2, Dopant::N_TYPE));
      contact.set_mobilities(70, 54);
      ModelOptions opts;
      opts["tau_n"] = "2e-9";
      opts["tau_p"] = "6e-10";
      RecombinationModelInterface* rm =
        RecombinationModelInterface::create("SRH", opts);
      contact.add_recombination_model(rm);
    }

    SemiconductorModel contact2;
    contact2.set_data_file("Si.dat");
    contact2.read_database(d);
    {
      contact2.add_dopant(new Dopant(5e19, 0.025, 2, Dopant::N_TYPE));
      contact2.set_mobilities(70, 54);
      ModelOptions opts;
      opts["tau_n"] = "2e-9";
      opts["tau_p"] = "6e-10";
      //opts.insert(pair<const string, string>("tau_n", "2e-9"));
      //opts.insert(pair<const string, string>("tau_p", "6e-10"));
      RecombinationModelInterface* rm =
        RecombinationModelInterface::create("SRH", opts);
      contact2.add_recombination_model(rm);
    }



    if (statistics == "FD")
    {
      sub.set_statistics(TiberCad::FERMIDIRAC);
      schottky.set_statistics(TiberCad::FERMIDIRAC);
      contact.set_statistics(TiberCad::FERMIDIRAC);
    }
    else
    {
      sub.set_statistics(TiberCad::BOLTZMANN);
      schottky.set_statistics(TiberCad::BOLTZMANN);
      contact.set_statistics(TiberCad::BOLTZMANN);
    }
    
    if (fully)
    {
      sub.set_coupling_type(BOTH);
      schottky.set_coupling_type(BOTH);
      contact.set_coupling_type(BOTH);
    }
    else
    {
      sub.set_coupling_type(ELECTRONS);
      schottky.set_coupling_type(ELECTRONS);
      contact.set_coupling_type(ELECTRONS);
    }


    ElementData element_data;
    {
      MeshData::const_elem_data_iterator it = meshdata.elem_data_begin();
      const MeshData::const_elem_data_iterator end = meshdata.elem_data_end();
      for ( ; it != end; ++it)
      {
        const Elem* elem = it->first;

        // every element needs to have a material assigned
        assert(meshdata.has_data(elem));

        int id = (int) meshdata(elem);

        switch (id)
        {
          case 2:
            element_data.set_data(elem, &schottky);
            break;
          case 3:
            element_data.set_data(elem, &contact);
            break;
          case 4:
            element_data.set_data(elem, &contact2);
            break;
          default:
            element_data.set_data(elem, &sub);
            break;
        }
      }
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
            set_boundary(boundary_data, nodes, &source, mesh);
            break;
          case 2:
            set_boundary(boundary_data, nodes, &gate, mesh);
            break;
          case 3:
            set_boundary(boundary_data, nodes, &drain, mesh);
            break;
        }
      }
    }


    DD::Device device(&mesh, &element_data, &boundary_data);
    bool device_integrity = device.check_integrity();
    if (device_integrity)
      cout << "Device ok.\n\n";
    else
      cout << "Device bad.\n\n";
    
  
    EquationSystems eqsys(mesh);
    DriftDiffusion dd(&device);
    dd.set_equation_systems(&eqsys);
    dd.init();

    DriftDiffusion::Options& params = dd.get_options();
    params.solver_params.nonlinear_max_iterations = 1;
    params.solver_params.linear_max_iterations = lin_max_it;
    params.solver_params.ls_maxstep = nonlin_ls_maxstep;
    params.solver_params.nonlinear_tolerance = nonlin_rtol;
    params.solver_params.ls_maxstep = nonlin_ls_maxstep;
    params.solver_params.linear_tolerance = lin_rtol;
    params.integration_order =
      static_cast<libMeshEnums::Order>(integration_order);
    params.max_refinement_steps = refinement_steps;
    params.refine_fraction = refine_frac;
    params.coarsen_fraction = coarsen_frac;

    //params.local_scaling = true;

    // mesh drawn in um
    params.mesh_units = mesh_units;

    dd.enable_mesh_refinement();
    
    //print_boundary_data(boundary_data, mesh);

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

    cout << "Material properties:\n";
    sub.print_info();

    vector<double> densities;
    vector<string> names;
    dd.build_band_edges(densities, names);
    TecplotIO(dd.get_mesh()).write_nodal_data("output/eq_bands.plt",
        densities, names);

    dd.build_densities(densities, names);
    TecplotIO(dd.get_mesh()).write_nodal_data("output/eq_densities.plt",
        densities, names);



    params.solver_params.nonlinear_max_iterations = nonlin_max_it;

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
        DriftDiffusion::Options& params = dd.get_options();
        int bkp = params.coupling;
        params.coupling = POISSON;
        params.solver_params.nonlinear_max_iterations = 100;
        dd.set_simulation_voltage("gate", 0.0);
        dd.set_simulation_voltage("drain", 0.0);
        dd.set_electron_fermi_level(0.0);
        dd.set_hole_fermi_level(0.0);
        dd.set_electric_potential(0.0);
        dd.guess_equilibrium();
        dd.solve();
        dd.remember_current_solution();
        params.coupling = bkp;
        params.solver_params.nonlinear_max_iterations = nonlin_max_it;

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


void set_boundary(BoundaryData& data, const vector<unsigned int>& nodes,
    ElectricalContact* desc, const Mesh& mesh)
{
  vector<unsigned int>::const_iterator n_it;
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


