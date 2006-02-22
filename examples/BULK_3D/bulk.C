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
    DriftDiffusionProperties& sc_model);

void set_boundary(BoundaryData& data, const vector<unsigned int>& nodes,
    BoundaryDescriptor& desc, const Mesh& mesh);

void print_boundary_data(const BoundaryData& data, const Mesh& mesh);



int main (int argc, char** argv)
{
  libMesh::init(argc, argv);
  {

    GetPot input("bulk.in");

    const string meshfile = input("meshfile", "");

    double temperature = input("temperature", 300.0);

    double start_voltage = input("start_voltage", 0.0);
    double stop_voltage = input("stop_voltage", 0.0);
    unsigned int voltage_steps = input("voltage_steps", 1);


    const string method = input("simulation_method", "NEWTON");
    const string statistics = input("statistics", "B");

    double nonlin_rtol = input("nonlinear_tolerance", 1e-9);
    double lin_rtol = input("linear_tolerance", 1e-12);
    int integration_order = input("integration_order", 5);
    int nonlin_max_it = input("nonlinear_max_it", 15);
    int lin_max_it = input("linear_max_it", 500);
    double nonlin_ls_maxstep =
      input("nonlinear_ls_maxstep", 0.025);

    double n_doping = input("n_doping", 1e15);
    double p_doping = input("p_doping", 0.0);
    double polarization_x = input("polarization_x", 0.0);
    double polarization_y = input("polarization_y", 0.0);
    
    double mesh_units = input("mesh_units", 0.0);


    unsigned int refinement_steps =
      input("max_refinement_steps", 0);
    double refine_frac = input("refine_fraction", 0.7);
    double coarsen_frac = input("coarsen_fraction", 0.3);

    vector<unsigned int> phys_reg_ID(1);
    phys_reg_ID[0] = 1; // p
    //phys_reg_ID[1] = 2; // n
    
    vector<unsigned int> BC_reg_ID(2);
    BC_reg_ID[0] = 3; // anode
    BC_reg_ID[1] = 4; // cathode

    
    int dim = 3;
    Mesh mesh(dim);
    MeshData_elements meshdata(mesh);
    meshdata.enable_compatibility_mode();
    
    cerr << "Read meshfile: " << meshfile << "\n";
    Read_MSH readmesh(meshfile, phys_reg_ID, BC_reg_ID, dim, mesh, meshdata);
    map<unsigned int, vector<unsigned int> > boundary_nodes;
    readmesh.get_BC_data(boundary_nodes);

    mesh.print_info();

    SimpleSemiconductorModel nside;    
    if (statistics == "FD")
      nside.set_statistics(TiberCad::FERMIDIRAC);
    else
      nside.set_statistics(TiberCad::BOLTZMANN);

    nside.add_recombination_model(SRH);

    nside.set_relative_permittivity(11.7);
    nside.set_valence_band_properties(-0.5, 0.81, 200);
    nside.set_conduction_band_properties(0.62, 1.18, 800);
    nside.set_SRH_parameters(1e-8, 1e-8);

    SimpleSemiconductorModel pside(nside);
    pside.set_SRH_parameters(1e-8, 1e-8);

    nside.set_n_dopant(Dopant(n_doping, 0.025, 2));
    pside.set_n_dopant(Dopant(n_doping, 0.025, 2));
    nside.set_p_dopant(Dopant(p_doping, 0.01, 4));
    pside.set_p_dopant(Dopant(p_doping, 0.01, 4));

    nside.calculate_equilibrium_properties(BOTH);
    pside.calculate_equilibrium_properties(BOTH);


    ElementData element_data;

    cerr << "setup element data\n";
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
          default:
            element_data.set_data(elem, &nside);
            break;
        }
      }
    }

    BoundaryDescriptor anode("anode");
    BoundaryDescriptor cathode("cathode");

    setup_boundary_desc(anode, nside);
    cerr << "p side:\n";
    cerr << pside.get_equilibrium_fermi_level() << " eV, " << "ni = " <<
      pside.get_intrinsic_density() << " n0 = " <<
      pside.get_equilibrium_electron_density() << " p0 = " << 
      pside.get_equilibrium_hole_density() << "\n";

    setup_boundary_desc(cathode, nside);
    cerr << "n side:\n";
    cerr << nside.get_equilibrium_fermi_level() << " eV, " << "ni = " <<
      nside.get_intrinsic_density() << " n0 = " <<
      nside.get_equilibrium_electron_density() << " p0 = " << 
      nside.get_equilibrium_hole_density() << "\n";



    BoundaryData boundary_data;
    {
      map<unsigned int, vector<unsigned int> >::const_iterator it =
        boundary_nodes.begin();
      const map<unsigned int, vector<unsigned int> >::const_iterator end =
        boundary_nodes.end();

      for ( ; it != end; ++it)
      {
        const vector<unsigned int>& nodes = it->second;

        if (it->first == BC_reg_ID[0])
          set_boundary(boundary_data, nodes, anode, mesh);
        else
          set_boundary(boundary_data, nodes, cathode, mesh);
      }
    }


    //print_boundary_data(boundary_data, mesh);


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

    //SemiconductorModel sc_model;
    //dd.set_semiconductor_model(&sc_model);

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

    if (method == "GUMMEL")
    {
      params.solver_method = DriftDiffusion::GUMMEL;
      params.max_gummel_iterations = 2;
    }
    else
      params.solver_method = DriftDiffusion::NEWTON;

    params.coupling = FULLYCOUPLED;

    params.approximation_order = FIRST;

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
      GMVIO(dd.get_mesh()).write_nodal_data("output/nodal_eq.gmv",
          dd.get_solution(), dd.get_variable_names());
      dd.build_densities(densities, names);
      GMVIO(dd.get_mesh()).write_nodal_data("output/densities_eq.gmv",
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
      dd.solve();
      cout << "done (nr. iterations: " << dd.get_n_nonlinear_iterations() <<
        ", final residual: " << dd.get_final_residual() << ")\n" << flush;
      ostringstream filename;
      filename << "output/nodal_" << *it << ".gmv";
      GMVIO(dd.get_mesh()).write_nodal_data(filename.str(),
          dd.get_solution(), dd.get_variable_names());
      ostringstream filename_d;
      filename_d << "output/densities_" << *it << ".gmv";
      vector<double> densities;
      vector<string> names;
      dd.build_densities(densities, names);
      GMVIO(dd.get_mesh()).write_nodal_data(filename_d.str(),
          densities, names);

      const map<const BoundaryDescriptor*, double>& curr =
        dd.get_boundary_currents();
      iv_char[*it] = vector<double>(3);
      iv_char[*it][0] = (*curr.find(&cathode)).second;
      iv_char[*it][1] = (*curr.find(&anode)).second;
      iv_char[*it][2] = dd.get_artificial_boundary_current();
      cerr << "    I = " << (iv_char[*it][1]) << " A\n";
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
        dd.solve(restart);
        cout << "done (nr. iterations: " << dd.get_n_nonlinear_iterations() <<
          ", final residual: " << dd.get_final_residual() << ")\n" << flush;
        ostringstream filename;
        filename << "output/nodal_" << *it << ".gmv";
        GMVIO(dd.get_mesh()).write_nodal_data(filename.str(),
            dd.get_solution(), dd.get_variable_names());
        ostringstream filename_d;
        filename_d << "output/densities_" << *it << ".gmv";
        vector<double> densities;
        vector<string> names;
        dd.build_densities(densities, names);
        GMVIO(dd.get_mesh()).write_nodal_data(filename_d.str(),
            densities, names);

        restart = false;
        const map<const BoundaryDescriptor*, double>& curr =
          dd.get_boundary_currents();
        iv_char[*it] = vector<double>(3);
        iv_char[*it][0] = (*curr.find(&cathode)).second;
        iv_char[*it][1] = (*curr.find(&anode)).second;
        iv_char[*it][2] = dd.get_artificial_boundary_current();
        cerr << "    I = " << (iv_char[*it][1]) << " A\n";
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

void setup_boundary_desc(BoundaryDescriptor& desc,
    DriftDiffusionProperties& sc)
{
  std::vector<double> coeff(3, 0);
  coeff[0] = 1.0;

  desc.set_coefficients("fermi_e", coeff);
  desc.set_coefficients("fermi_h", coeff);

  coeff[2] = sc.get_equilibrium_fermi_level();

  cout << desc.get_id() << ": " << coeff[2] << "\n";
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
      cout << "  (" << n(0) << ", " << n(1) << ", " << n(2) << ")";
    }
    cout << "\n";
    /*
    cout << "  boundary values: ";
    cout << "a = " << val[0] << " ";
    cout << "b = " << val[1] << " ";
    cout << "c = " << val[2] << " ";
    cout << "\n\n";
    */

    ++it;
  }
}
