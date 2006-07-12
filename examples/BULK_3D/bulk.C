// the following _HAS_ to be included first
#include "Read_MSH.h"

#include "ElementData.h"
#include "OhmicContact.h"
#include "FermiLevelPinning.h"
#include "SchottkyContact.h"
#include "BoundaryData.h"
#include "DDevice.h"
#include "DriftDiffusion.h"
#include "SemiconductorModel.h"
#include "RecombinationModelInterface.h"

#include "mesh.h"
#include "equation_systems.h"
#include "mesh_modification.h"
#include "mesh_data.h"
#include "mesh_generation.h"
#include "elem.h"
#include "getpot.h"
#include "GMVIO_cell.h"

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

    GetPot input("bulk.in");

    const string meshfile = input("meshfile", "");
    const int dim = input("dimension", 3);

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
    phys_reg_ID[0] = 1; 
    
    vector<unsigned int> BC_reg_ID(2);
    BC_reg_ID[0] = 3; // anode
    BC_reg_ID[1] = 4; // cathode

    
    Mesh mesh(dim);
    MeshData_elements meshdata(mesh);
    meshdata.enable_compatibility_mode();
    
    cerr << "Read meshfile: " << meshfile << "\n";
    Read_MSH readmesh(meshfile, phys_reg_ID, BC_reg_ID, dim, mesh, meshdata);
    map<unsigned int, vector<unsigned int> > boundary_nodes;
    readmesh.get_BC_data(boundary_nodes);

    mesh.print_info();

    Dummy d;
    SemiconductorModel nside;    
    nside.set_data_file("Si.dat");
    nside.read_database(d);

    if (statistics == "FD")
      nside.set_statistics(TiberCad::FERMIDIRAC);
    else
      nside.set_statistics(TiberCad::BOLTZMANN);

    ModelOptions opts;
    opts["tau_n"] = "1e-7";
    opts["tau_p"] = "3e-8";
    RecombinationModelInterface* recomb_model =
      RecombinationModelInterface::create("SRH", opts);
    nside.add_recombination_model(recomb_model);


    nside.set_mobilities(800, 200);

    nside.add_dopant(new Dopant(n_doping, 0.025, 2, Dopant::N_TYPE));
    nside.add_dopant(new Dopant(p_doping, 0.01, 4, Dopant::P_TYPE));


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

    OhmicContact anode("anode");
    //SchottkyContact anode("anode");
    //anode.set_schottky_barrier(0.8);
    //FermiLevelPinning anode("anode");
    //anode.set_pinning(1.0);
    OhmicContact cathode("cathode");
    //FermiLevelPinning cathode("cathode");
    //cathode.set_pinning(1.0);

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
      dd.remember_current_solution();
      cout << "done (nr. iterations: " << dd.get_n_nonlinear_iterations() <<
        ", final residual: " << dd.get_final_residual() << ")\n" << flush;
      dd.build_band_edges(densities, names);
      GMVIO(dd.get_mesh()).write_nodal_data("output/bands_eq.gmv",
          densities, names);
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

    // calculate for electron and holes
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
      //if (*it >= 1.7)
      //  dd.enable_mesh_refinement();
      
      dd.solve();
      
      cout << "done (nr. iterations: " << dd.get_n_nonlinear_iterations() <<
        ", final residual: " << dd.get_final_residual() << ")\n" << flush;
      ostringstream filename;
      vector<double> densities;
      vector<string> names;
      ostringstream filename_b;
      filename_b << "output/bands_" << *it << ".gmv";
      dd.build_band_edges(densities, names);
      GMVIO(dd.get_mesh()).write_nodal_data(filename_b.str(),
          densities, names);
      ostringstream filename_d;
      filename_d << "output/densities_" << *it << ".gmv";
      dd.build_densities(densities, names);
      GMVIO(dd.get_mesh()).write_nodal_data(filename_d.str(),
          densities, names);
      ostringstream filename_f;
      filename_f << "output/field_" << *it << ".gmv";
      dd.build_electric_field(densities, names);
      GMVIO_cell(dd.get_mesh()).write_ascii_cell_data(filename_f.str(),
          densities, names);
      ostringstream filename_c;
      filename_c << "output/current_" << *it << ".gmv";
      dd.build_current_density(densities, names);
      GMVIO_cell(dd.get_mesh()).write_ascii_cell_data(filename_c.str(),
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
      dd.set_to_remembered_solution();
      do
      {
        --it;
        dd.set_simulation_voltage("anode", *it);
        cout << " Solving U = " << *it << " V ...\n" << flush;
        dd.enable_mesh_refinement();
        dd.solve();
        cout << "done (nr. iterations: " << dd.get_n_nonlinear_iterations() <<
          ", final residual: " << dd.get_final_residual() << ")\n" << flush;
        vector<double> densities;
        vector<string> names;
        ostringstream filename_b;
        filename_b << "output/bands_" << *it << ".gmv";
        dd.build_band_edges(densities, names);
        GMVIO(dd.get_mesh()).write_nodal_data(filename_b.str(),
            densities, names);
        ostringstream filename_d;
        filename_d << "output/densities_" << *it << ".gmv";
        dd.build_densities(densities, names);
        GMVIO(dd.get_mesh()).write_nodal_data(filename_d.str(),
            densities, names);
        ostringstream filename_f;
        filename_f << "output/field_" << *it << ".gmv";
        dd.build_electric_field(densities, names);
        GMVIO_cell(dd.get_mesh()).write_ascii_cell_data(filename_f.str(),
            densities, names);
        ostringstream filename_c;
        filename_c << "output/current_" << *it << ".gmv";
        dd.build_current_density(densities, names);
        GMVIO_cell(dd.get_mesh()).write_ascii_cell_data(filename_c.str(),
            densities, names);

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
    ElectricalContact *desc, const Mesh& mesh)
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


