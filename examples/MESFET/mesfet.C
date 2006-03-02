// the following _HAS_ to be included first
#include "Read_MSH.h"

#include "ElementData.h"
#include "BoundaryDescriptor.h"
#include "BoundaryData.h"
#include "DDevice.h"
#include "DriftDiffusion.h"
#include "SimpleSemiconductorModel.h"

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

void setup_boundary_desc(BoundaryDescriptor& desc,
    const DriftDiffusionProperties& sc_model);

void set_boundary(BoundaryData& data, const vector<unsigned int>& nodes,
    BoundaryDescriptor& desc, const Mesh& mesh);

void print_boundary_data(const BoundaryData& data, const Mesh& mesh);

void sweep_drain(double stop, int steps, DriftDiffusion& dd, double vg,
    bool restart = false);



int main (int argc, char** argv)
{
  libMesh::init(argc, argv);
  {

    GetPot input("mesfet.in");

    string meshfile = input("meshfile", "");
    const double mesh_units = input("mesh_units", 1e-4);

    int fully = input("fully_coupled", 1);

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

    vector<unsigned int> phys_reg_ID(3);
    phys_reg_ID[0] = 101; // bulk
    phys_reg_ID[1] = 102; // contacts
    phys_reg_ID[2] = 103; // doped
    
    vector<unsigned int> BC_reg_ID(3);
    BC_reg_ID[0] = 4; // source
    BC_reg_ID[1] = 5; // gate
    BC_reg_ID[2] = 6; // drain

    int dim = 2;
    Mesh mesh(dim);
    MeshData_elements meshdata(mesh);
    meshdata.enable_compatibility_mode();
    
    Read_MSH readmesh(meshfile, phys_reg_ID, BC_reg_ID, dim, mesh, meshdata);
    map<unsigned int, vector<unsigned int> > boundary_nodes;
    readmesh.get_BC_data(boundary_nodes);

    mesh.print_info();

    SimpleSemiconductorModel sub;    

    if (statistics == "FD")
      sub.set_statistics(TiberCad::FERMIDIRAC);
    else
      sub.set_statistics(TiberCad::BOLTZMANN);

    sub.add_recombination_model(SRH);

    sub.set_relative_permittivity(11.7);
    sub.set_valence_band_properties(-0.5, 0.81, 200);
    sub.set_conduction_band_properties(0.62, 1.18, 800);
    sub.set_SRH_parameters(1e-7, 1e-7);

    SimpleSemiconductorModel schottky(sub);
    
    schottky.set_n_dopant(Dopant(1e18, 0.025, 2));
    schottky.set_SRH_parameters(1e-8, 1e-8);

    SimpleSemiconductorModel contact(sub);

    contact.set_n_dopant(Dopant(5e19, 0.025, 2));
    contact.set_SRH_parameters(1e-9, 1e-9);

    const Elem* elem = meshdata.elem_data_begin()->first;
    sub.reinit(elem);
    schottky.reinit(elem);
    contact.reinit(elem);
    sub.calculate_equilibrium_properties();
    schottky.calculate_equilibrium_properties();
    contact.calculate_equilibrium_properties();

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
          case 102:
            element_data.set_data(elem, &contact);
            break;
          case 103:
            element_data.set_data(elem, &schottky);
            break;
          default:
            element_data.set_data(elem, &sub);
            break;
        }
      }
    }


    BoundaryDescriptor source("source");
    BoundaryDescriptor gate("gate");
    BoundaryDescriptor drain("drain");
    setup_boundary_desc(source, contact);
    setup_boundary_desc(drain, contact);
    cerr << "contact:\n";
    cerr << contact.get_equilibrium_fermi_level() << " eV, " << "ni = " <<
      contact.get_intrinsic_density() << " n0 = " <<
      contact.get_equilibrium_electron_density() << " p0 = " << 
      contact.get_equilibrium_hole_density() << "\n";
    setup_boundary_desc(gate, schottky);
    cerr << "schottky:\n";
    cerr << schottky.get_equilibrium_fermi_level() << " eV, " << "ni = " <<
      schottky.get_intrinsic_density() << " n0 = " <<
      schottky.get_equilibrium_electron_density() << " p0 = " << 
      schottky.get_equilibrium_hole_density() << "\n";

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
          case 4:
            set_boundary(boundary_data, nodes, source, mesh);
            break;
          case 5:
            set_boundary(boundary_data, nodes, gate, mesh);
            break;
          case 6:
            set_boundary(boundary_data, nodes, drain, mesh);
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
    
  
    DriftDiffusion dd(&device);

    DriftDiffusion::Options& params = dd.get_options();
    params.solver_params.nonlinear_max_iterations = nonlin_max_it;
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

    // mesh drawn in um
    params.mesh_units = mesh_units;

    dd.enable_mesh_refinement();
    
    //print_boundary_data(boundary_data, mesh);

    dd.set_simulation_voltage("source", 0.0);
    dd.set_simulation_voltage("gate", 0.0);
    dd.set_simulation_voltage("drain", 0.0);
    
    if (characteristic)
    {
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
      GMVIO(dd.get_mesh()).write_nodal_data("output/eq_potentials.gmv",
          dd.get_solution(), dd.get_variable_names());
      vector<double> densities;
      vector<string> names;
      dd.build_densities(densities, names);
      GMVIO(dd.get_mesh()).write_nodal_data("output/eq_densities.gmv",
          densities, names);

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
        const map<const BoundaryDescriptor*, double>& curr =
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
          const map<const BoundaryDescriptor*, double>& curr =
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


/*
   const map<const BoundaryDescriptor*, double>& curr =
   dd.get_boundary_currents();
   double is = curr.find(&source)->second;
   double ig = curr.find(&gate)->second;
   double id = curr.find(&drain)->second;
   double check = dd.get_artificial_boundary_current();

   cout << "Currents:\n";
   cout << " is = " << is << "\n";
   cout << " ig = " << ig << "\n";
    cout << " id = " << id << "\n";
    cout << " (is + ig + id = " << is + ig + id << ")\n";
    cout << " (artificial boundary: i = " << check << ")\n";
*/


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
    //dd.enable_mesh_refinement();
    dd.solve(restart);
    if (remember)
    {
      dd.remember_current_solution();
    }
    remember = false;
    restart = false;
    const map<const BoundaryDescriptor*, double>& curr =
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
    GMVIO(dd.get_mesh()).write_nodal_data("output/potentials"+f.str(),
        dd.get_solution(), dd.get_variable_names());
    dd.build_densities(densities, names);
    GMVIO(dd.get_mesh()).write_nodal_data("output/densities"+f.str(),
        densities, names);
  }
  
  file.close();
}


void setup_boundary_desc(BoundaryDescriptor& desc,
    const DriftDiffusionProperties& sc)
{
  std::vector<double> coeff(3, 0);
  coeff[0] = 1.0;

  desc.set_coefficients("fermi_e", coeff);
  desc.set_coefficients("fermi_h", coeff);

  if (desc.get_id() != "gate")
  {
    coeff[2] = sc.get_equilibrium_fermi_level();
  }
  else
  {
    const SimpleSemiconductorModel& m =
      static_cast<const SimpleSemiconductorModel&>(sc);
    double ec = m.get_conduction_band_edge();
    double barrier = 0.8;
    coeff[2] = ec - barrier;    
  }

  cout << coeff[2] << "\n";
  desc.set_coefficients("potential", coeff);

}


void set_boundary(BoundaryData& data, const vector<unsigned int>& nodes,
    BoundaryDescriptor& desc, const Mesh& mesh)
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
          if (find(n_begin, n_end, side->node(i) + 1) == n_end)
            found = false;
        }
        if (found)
          data.set_data(BoundaryData::ElementSide(elem, s), &desc);
      }
    }
  }
}


void print_boundary_data(const BoundaryData& data, const Mesh& mesh)
{

  BoundaryData::const_iterator it = data.sides_begin();
  const BoundaryData::const_iterator end = data.sides_end();
  while (it != end)
  {
    const BoundaryData::ElementSide& s = it->first;
    const BoundaryDescriptor* desc = it->second;
    const vector<double>& val = *(desc->get_coefficients("potential"));

    cout << desc->get_id();
    cout << " nodes:";
    AutoPtr<Elem> side = s.first->build_side(s.second);
    for (int i = 0; i < side->n_nodes(); i++)
    {
      Point& n = side->point(i);
      cout << "  (" << n(0) << ", " << n(1) << ")";
    }
    cout << "\n";
    cout << "  boundary values: ";
    cout << "a = " << val[0] << " ";
    cout << "b = " << val[1] << " ";
    cout << "c = " << val[2] << " ";
    cout << "\n\n";

    ++it;
  }
}
