// $Id$

#include "EigenvalueProblem.h"
#include "SimulationEnvironment.h"
#include "AtomisticStructure.h"
#include "BulkCrystal.h"
#include "Constants.h"
#include "Messages.h"
#include "DataOutput.h"

#include "SimulationOptions.h"
#include <petsc_matrix.h>


#include <petsc_matrix.h>

#include "libmesh/elem.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/mesh_tools.h"

#include <boost/shared_ptr.hpp>

using namespace std;
using namespace libMesh;


namespace
{
  bool compare_points(const Point& a, const Point& b)
  {
    return(a.absolute_fuzzy_equals(b, 1e-5));
  }
}

void
EigenvalueProblem::initialize_solution_container(size_t num_solutions)
{
  _solution.clear();
  _solution.resize(num_solutions);
}


void EigenvalueProblem::init_kspace(const ModelOptions& opt)
{
   ModelOptions kopts = parse_kspace_options(opt);

   _kspace = new Kspace(kopts, get_communicator());

   if(_kspace==NULL)
     throw InitFailedException("Could not initialize k-space");
   else
     Messages::info("k-space initialized");

}
 
ModelOptions EigenvalueProblem::parse_kspace_options(const ModelOptions& opts)
{
  ModelOptions kopts(opts);

  //kopts.set_option("mesh_units", get_mesh_units());
  unsigned int mesh_dim = get_mesh().mesh_dimension();
  
  bool x_periodic = get_option("x-periodicity", false);
  bool y_periodic = get_option("y-periodicity", (mesh_dim < 2));
  bool z_periodic = get_option("z-periodicity", (mesh_dim < 3));

  if (opts.find_option("k_path") || opts.find_option("k-path"))
  {	  
    std::string kpath = opts.get_option("k-path","");
    kpath = opts.get_option("k_path",kpath);    
    kopts.set_option("k-path",kpath);
    //std::vector<unsigned int>  num_nodes(1,20);
    //kopts.set_option("number_of_nodes",num_nodes);
    //ModelOptions newopts;
    //newopts.set_option("output_format","grace");
    //set_options( newopts );
  }

  // these are the real space lattice vectors, in nm
  // why pi? Because then the default max k becomes 1 ( = 2*pi/(2*a) ), and
  // k_max can be interpreted in nm^-1
  // we make it slightly anisotropic so it will be detected as orthorhombic
  RealVectorValue a(M_PI, 0.0, 0.0), b(0.0, 1.01*M_PI, 0.0), c(0.0, 0.0, 1.005*M_PI);
  auto bbox = get_environment().get_bounding_box();
  if (get_option("x-periodicity", false) && (mesh_dim > 0))
    a(0) = (bbox.second(0) - bbox.first(0)) * get_mesh_units() * 1e9;
  if (get_option("y-periodicity", false) && (mesh_dim > 1))
    b(1) = (bbox.second(1) - bbox.first(1)) * get_mesh_units() * 1e9;
  if (get_option("z-periodicity", false) && (mesh_dim > 2))
    c(2) = (bbox.second(2) - bbox.first(2)) * get_mesh_units() * 1e9;


  // if there is an atomistic structure, we can take the lattice vectors from it
  // (they come in Angstrom!)
  if (get_atomistic_structure() != NULL)
  {
    get_atomistic_structure()->get_lattice_vectors(a, b, c);
    a *= 0.1;
    b *= 0.1;
    c *= 0.1;

    x_periodic = get_atomistic_structure()->is_periodic(0);
    y_periodic = get_atomistic_structure()->is_periodic(1);
    z_periodic = get_atomistic_structure()->is_periodic(2);
  }

  // find k-space dimension, considering real space periodicities
  unsigned int k_dim = 0;
  if (x_periodic)
    k_dim++;
  if (y_periodic)
    k_dim++;
  if (z_periodic)
    k_dim++;

  //k_dim = min(k_dim, 3u);

  if (opts.find_option("k_space_dimension"))
    k_dim = opts.get_option("k_space_dimension", k_dim);

  kopts.set_option("k_space_dimension", k_dim);

  // we skip this if dimension are provided explicitly
  if (!opts.find_option("r1"))
  {

  // for now, we have the old approach for continuum models,
  // while we use directly the periodicity info for atomistic models
    switch (k_dim)
    {
      case 1:
        if (get_atomistic_structure() != NULL)
        {
          if (x_periodic)
            c = a;
          else if (y_periodic)
            c = b;
        }
        else if (mesh_dim == 3)
        {
          if (x_periodic)
            c = a;
          else if (y_periodic)
            c = b;
        }
        kopts.set_option("r1", c);
        break;

      case 2:
        if (get_atomistic_structure() != NULL)
        {
          if (!y_periodic)
            b = a;
          else if (!z_periodic)
          {
            c = b;
            b = a;
          }

        }
        else if (mesh_dim == 3)
        {
          if (!z_periodic)
          {
            c = b;
            b = a;
          }
          else if (x_periodic)
          {
            b = a;
          }
        }
        else if (mesh_dim == 2)
        {
          if (x_periodic)
            b = a;
        }
        kopts.set_option("r1", b);
        kopts.set_option("r2", c);
        break;

      case 3:
        kopts.set_option("r1", a);
        kopts.set_option("r2", b);
        kopts.set_option("r3", c);
        break;

      default:
        break;
    }
  }


  kopts.set_option("mesh_order", opts.get_option("mesh_order", "first"));

  return kopts;
}





Point EigenvalueProblem::get_k_point(bool relative_coord) const
{
  Point kp(_k_vector[0], _k_vector[1], _k_vector[2]);
  if (relative_coord)
  {
    _kspace->inverse_transform(kp);
  }

  return kp;
}




void EigenvalueProblem::do_calculate_density_at_k(DofField&)
{
}

void EigenvalueProblem::do_assemble(const ModelOptions&)
{
}


void EigenvalueProblem::get_H_csr(std::vector<libMesh::Complex>&, std::vector<int>&,
        std::vector<int>&) const
{
}


void EigenvalueProblem::get_S_csr(std::vector<libMesh::Complex>&, std::vector<int>&,
        std::vector<int>&) const
{
}

void EigenvalueProblem::print_H(const std::string&) const
{
}


double EigenvalueProblem::get_band_edge(const std::string&)
{
  return 0;
}


void EigenvalueProblem::compute_dispersion(const ModelOptions& opts)
{
  ModelOptions kopts(opts);
  string refmat_s = kopts.get_option("unfold_to", "");
  kopts.delete_option("unfold_to");

  Material* refmat = nullptr;
  if (!refmat_s.empty())
  {
    std::set<ID> ids;
    get_environment().get_device().get_active_region_ids(refmat_s, ids);
    if (ids.size() != 0)
    {
      refmat = get_environment().get_device().get_material(*ids.begin());
    }
  }
  else
  {
    ModelOptions::const_submodel_iterator it =
        kopts.submodels_begin("unfold_to");
    if (it != kopts.submodels_end("unfold_to"))
    {
      const ModelOptions& refmat_opts = it->second;
      refmat = Material::create(refmat_opts.get_name(), refmat_opts);
      refmat->init();
    }
  }

  // this will contain the necessary k-points
  vector<Point> Kpoints;
  // these map k points from SC to PC, in case of unfolding
  vector<set<unsigned int>> k_to_K;
  vector<set<unsigned int>> K_to_k;

  const MeshBase* kmesh = nullptr;
  if (refmat != nullptr)
  {
    Messages::newline();
    Messages::info("Dispersion will be unfolded to BZ of " +
        refmat->get_name());

    ModelOptions scopts(kopts);
    scopts.delete_option("k-path");
    scopts.delete_option("k_path");
    //scopts.set_option("write_k_mesh", "SC");

    Messages::newline();
    Messages::info("Initialize super cell k-space");
    init_kspace(scopts);

    // backup the SC k-space
    Kspace* sc_kspace = _kspace;

    Messages::newline();
    Messages::info("Initialize primitive cell k-space");

    // build the BulkCrystal object
    BulkCrystal* bulk = BulkCrystal::create(refmat);
    bulk->init();

    kopts.set_option("k_space_dimension", 3);

    RealVectorValue a, b, c;
    bulk->get_lattice_vectors(a, b, c);
    a *= 0.1;
    b *= 0.1;
    c *= 0.1;

    kopts.set_option("r1", a);
    kopts.set_option("r2", b);
    kopts.set_option("r3", c);


    // _kspace contains now the k-path in the PC k-space
    init_kspace(kopts);

    // this is the reduced BZ of the PC k-space
    kopts.delete_option("k-path");
    kopts.delete_option("k_path");
    //kopts.set_option("write_k_mesh", "PC");
    kopts = parse_kspace_options(kopts);
    Kspace pc_kspace(kopts, get_communicator());
    const MeshBase* pc_mesh = pc_kspace.get_k_mesh();

    // get the bounding box, to obtain BZ radius
    MeshTools::BoundingBox pc_bb =
        MeshTools::bounding_box(*pc_mesh);
    pc_bb.max() += Point(1e-2, 1e-2, 1e-2);
    pc_bb.min() -= Point(1e-2, 1e-2, 1e-2);
    double pc_diam = 2 * pc_bb.max().norm();

    // for debugging
    //ofstream of("test.dat");

    typedef vector<int> LPoint;

    // the lattice points in the SC lattice for
    // the unfolding vectors
    vector<LPoint> G;

    pc_kspace.get_basis(a, b, c);
    //cerr << "PC : " << endl;
    //cerr << a << endl;
    //cerr << b << endl;
    //cerr << c << endl;

    // the SC reciprocal basis
    sc_kspace->get_basis(a, b, c);
    //cerr << "SC : " << endl;
    //cerr << a << endl;
    //cerr << b << endl;
    //cerr << c << endl;


    // these are absolute upper limits
    int Na = ceil(0.5 * pc_diam / a.norm() + 1);
    int Nb = ceil(0.5 * pc_diam / b.norm() + 1);
    int Nc = ceil(0.5 * pc_diam / c.norm() + 1);

    for (int i = -Na; i < Na; ++i)
    {
      for (int j = -Nb; j < Nb; ++j)
      {
        for (int l = -Nc; l < Nc; ++l)
        {

          LPoint g = {i, j, l};

          Point gg(g[0], g[1], g[2]);

          sc_kspace->transform_point(gg);

          if (gg.norm() < (pc_diam + 1e-5))
          {
            //if (gg(0) < 1e-6)
            //  of << gg(1) << " " << gg(2) << endl;

            // to store equivalent points;
            vector<Point> points;
            pc_kspace.equivalent_points(gg, points, false);

            for (auto&& p : points)
            {
              if (pc_bb.contains_point(p))
              {
                bool found = false;
                for (unsigned int el = 0; el < pc_mesh->n_elem(); ++el)
                {
                  const Elem* elem = pc_mesh->elem_ptr(el);
                  if (elem->close_to_point(p, 1e-3))
                  {
                    G.push_back(g);

                    found = true;
                    break;
                  }
                }

                if (found)
                  break;
              }
            }
          }
        }
      }
    }


    // sort out repeated points, checking with the translation group
    {
      const vector<Point> star =
          { Point(1, 0, 0), Point(0, 1, 0), Point(0, 0, 1),
            Point(1, 1, 0), Point(0, 1, 1), Point(1, 0, 1),
            Point(1, -1, 0), Point(0, 1, -1), Point(1, 0, -1),
            Point(-1, 1, 1), Point(1, -1, 1), Point(1, 1, -1),
            Point(1, 1, 1),
          };

      set<int> ids;
      for (unsigned int i = 0; i < G.size(); ++i)
      {
        if (!ids.count(i))
        {
          Point g(G[i][0], G[i][1], G[i][2]);
          sc_kspace->transform_point(g);
          Point p(g);
          pc_kspace.inverse_transform(p);

          for (unsigned int j = 0; j < G.size(); ++j)
          {
            if ((i != j) && !ids.count(j))
            {
              Point g2(G[j][0], G[j][1], G[j][2]);
              sc_kspace->transform_point(g2);
              Point p2(g2);
              pc_kspace.inverse_transform(p2);

              bool equal = p.absolute_fuzzy_equals(p2, 1e-3);

              for (unsigned int s = 0; (s < star.size()) && !equal; ++s)
                equal |= p.absolute_fuzzy_equals(p2 + star[s], 1e-3);

              if (equal)
              {
                ids.insert(j);
              }
            }
          }
        }
      }

      vector<LPoint> Gtmp;
      for (unsigned int i = 0; i < G.size(); ++i)
      {
        if (!ids.count(i))
          Gtmp.push_back(G[i]);
      }

      G = Gtmp;
    }

    // for debugging
    //ofstream of2("test2.dat");


    ostringstream os;
    os << "Found " << G.size() << " G vectors unfolding the supercell "
        << "onto the primitive cell";
    Messages::info(os.str());

    /*
    for (int i = 0; i < G.size(); ++i)
    {
      Point gg(G[i][0], G[i][1], G[i][2]);
      sc_kspace->transform_point(gg);
      Point p(gg);
      pc_kspace.inverse_transform(p);
      cerr << "  (" << G[i][0] << ", " << G[i][1] << ", " << G[i][2] << ") " << p << endl;
      of2 << gg(1) << " " << gg(2) << " " << gg(0) << endl;
    }
    cerr << endl;
    */

    kmesh = _kspace->get_k_mesh();
    unsigned int num_k_points = kmesh->n_nodes();
    k_to_K.resize(num_k_points);

    // now construct all K (considering duplicates),
    // and construct K->k and k->K tables
    for (unsigned int i = 0; i < num_k_points; i++)
    {
      const Point& k = kmesh->point(i);

      for (unsigned int j = 0; j < G.size(); ++j)
      {
        Point gg(G[j][0], G[j][1], G[j][2]);
        sc_kspace->transform_point(gg);
        Point K = k - gg;

        vector<Point> eq_K;
        sc_kspace->equivalent_points(K, eq_K);

        unsigned int foundK = 0;
        for (unsigned int eq = 0 ; eq < eq_K.size(); ++eq)
        {
          for (foundK = 0; foundK < Kpoints.size(); ++foundK)
          {
            if (Kpoints[foundK].absolute_fuzzy_equals(eq_K[eq], 1e-3))
            {
              break;
            }
          }

          if (foundK < Kpoints.size())
            break;
        }

        if (foundK == Kpoints.size())
        {
          Kpoints.push_back(eq_K[0]);
          K_to_k.resize(K_to_k.size() + 1);
        }

        K_to_k[foundK].insert(i);
        k_to_K[i].insert(foundK);

      }

    }

    for (unsigned int i = 0; i < K_to_k.size(); ++i)
    {
      cerr << i << " : ";
      for (auto&& j : K_to_k[i])
        cerr << j << " ";
      cerr << endl;
    }
    cerr << endl;

    for (unsigned int i = 0; i < k_to_K.size(); ++i)
    {
      cerr << i << " : ";
      for (auto&& j : k_to_K[i])
        cerr << j << " ";
      cerr << endl;
    }
    cerr << endl;



    os.str("");
    os << num_k_points << " primitive cell k points are folded into "
        << Kpoints.size() << " supercell k points";
    Messages::info(os.str());
    Messages::newline();

    /*
    ofstream of_k("kpoints.dat");
    cerr << endl << "K points: " << endl;
    for (int i = 0; i < Kpoints.size(); ++i)
    {
      Point gg(Kpoints[i]);
      of_k << gg(1) << " " << gg(2) << endl;
      sc_kspace->inverse_transform(gg);
      cerr << Kpoints[i] << "  " << gg << endl;
    }
    cerr << endl;

    for (int i = 0; i < k_to_K.size(); ++i)
    {
      cerr << kmesh->point(i) << endl;
      for (int j = 0; j < k_to_K[i].size(); ++j)
        cerr << "  " << k_to_K[i][j] << endl;
    }
    */



    // does crash, why?
    //delete sc_kspace;
  }
  else
  {
    init_kspace(kopts);
    kmesh = _kspace->get_k_mesh();

    unsigned int number_of_k_points = kmesh->n_nodes();
    Kpoints.reserve(number_of_k_points);
    K_to_k.resize(number_of_k_points);
    for (unsigned int i = 0; i < number_of_k_points; i++)
    {
      const Point&  k_point = kmesh->point(i);
      Kpoints.push_back(k_point);
      K_to_k[i].insert(i);
    }
  }


  // now calculate for the K points

  unsigned int number_of_k_points = Kpoints.size();
  unsigned int number_of_eigs;

  {
    unsigned int i = 0;
    const Point&  k_point = Kpoints[i];

    solve_for_kpoint(k_point);
    number_of_eigs = get_num_states();

    std::vector<double> temp(number_of_eigs);
    _dispersion.resize(number_of_k_points, temp);

    for (unsigned int j = 0 ; j <  _dispersion[0].size(); j++)
      _dispersion[0][j] = _solution[j].eigen_energy;
  }

  for (unsigned int i = 1; i < number_of_k_points; i++)
  {
    const Point&  k_point = Kpoints[i];

    solve_for_kpoint(k_point);
    number_of_eigs = get_num_states();

    for (unsigned int j = 0 ; j < _dispersion[i].size() ; j++)
      _dispersion[i][j] = _solution[j].eigen_energy;
  }


}



void
EigenvalueProblem::plot_dispersion(const std::string& filename)
{
 
  const MeshBase* kmesh = _kspace->get_k_mesh();
  unsigned int number_of_k_points = kmesh->n_nodes();

    std::vector<std::string> formats;
    get_output_format(formats);


    short kdim =  _kspace->dimension();

    if (_kspace->is_k_path())
    {
      formats.resize(1);
      formats[0] = "grace";
    }



    for(short k=0; k<formats.size();k++)
    {

      std::string format = formats[k];

      if (!(_kspace->is_k_path()) &&
          (format == "grace") && (kdim > 1))
      {
        format = "vtk";
      }

      std::vector<double> results;
      std::vector<std::string> names;

      unsigned int number_of_eigs = _dispersion[0].size();
      names.resize(number_of_eigs);

      unsigned int number_of_k_points = kmesh->n_nodes();
      results.resize( number_of_eigs * number_of_k_points );

      for (unsigned int i = 0; i < number_of_eigs ; i++)
      {
        std::ostringstream i_str;
        //The states are numbered starting from 0
        i_str << "state_number_" << i;
        names[i] = i_str.str();

        for (unsigned int j = 0; j < number_of_k_points ; j++)
          results[number_of_eigs * j + i] = _dispersion[j][i];
      }


      DataOutput data_output(*kmesh, format);
      data_output.set_output_directory(get_output_directory());
      //data_output.set_filename(filename);

      data_output.write_nodal_data(filename, results, names);

    }
}


void
EigenvalueProblem::plot_globaldata(void)
{

  string outdir = get_output_directory();

  string filename(outdir + "/" + get_output_filename() + ".dat");

  ofstream file;
  file.open(filename.c_str());

  if (file.good())
  {
    // header
    file << "# " << get_type() << " eigenstates (" << get_name() << ")\n";

    file << "# Index" << setw(9)<< "Particle" << setw(12) << "EigenEnergy"
         << setw(15) << "Occupation"
         << setw(12) << "FermiLevel" << setw(12) << "Temperature" << "\n";

    for (unsigned int i = 0; i < _solution.size(); i++)
    {
        file << setw(7) << i << setw(8) << _solution[i].particle
             << setw(14) << _solution[i].eigen_energy << " "
             << setw(14) << get_population(i) << " "
             << setw(14) << _solution[i].electro_chem_pot << " "
             << setw(11) << _solution[i].temperature << "\n";
    }
  }

}





void
EigenvalueProblem::process_element(const Elem* elem, unsigned int entryside,
    vector<vector<eigen_problem_solution>>& ordered_solutions)
{
  ofstream of("test.dat", ofstream::app);

  bool already_done = true;

  // choose reference node
  unsigned int ref_node;
  for (unsigned int n = 0; n < elem->n_nodes(); ++n)
  {
    if (elem->is_node_on_side(n, entryside))
    { 
      ref_node = elem->node(n);
      if (ordered_solutions[ref_node].empty())
      {
        const Point& k_point = elem->point(n);

        solve_for_kpoint(k_point);
        int number_of_eigs = get_num_states();
        ordered_solutions[ref_node].resize(number_of_eigs);

        for (unsigned int j = 0 ; j < number_of_eigs; j++)
          ordered_solutions[ref_node][j] = _solution[j];
      }

      break;
    }
  }

  // the number of the reference solutions
  // (should be usually the same as the solution size)
  int ref_size = ordered_solutions[ref_node].size();

  for (unsigned int n = 0; n < elem->n_nodes(); ++n)
  {
    unsigned int node_id = elem->node(n);
    if (ordered_solutions[node_id].empty())
    {
      const Point& k_point = elem->point(n);

      solve_for_kpoint(k_point);
      int number_of_eigs = get_num_states();
      // we take the maximum only to make everything crash if the two numbers
      // do not correspond. However, this is usually sign of a badly posed
      // simulation setup
      ordered_solutions[node_id].resize(max(number_of_eigs, ref_size));


      //of << node_id << " ";
      //k_point.write_unformatted(of, true);

      set<unsigned int> ids;

      for (unsigned int k = 0 ; k < ref_size; k++)
      {
        unsigned int idx = k;
        double max_sp = 0;
        //cerr << k << " : ";
        for (unsigned int j = 0 ; j < number_of_eigs; j++)
        {
          if (!ids.count(j))
          {
            double proj =
                abs(scalar_product(ordered_solutions[ref_node][k], _solution[j]));
            //cerr << j << " - " << proj << " ";
            if (proj > max_sp)
            {
              max_sp = proj;
              idx = j;
            }
          }
        }
        //cerr << endl;
        ids.insert(idx);
        //cerr << idx << " " << max_sp << endl;
        ordered_solutions[node_id][k] = _solution[idx];
      }

      already_done = false;
    }
  }

  // go into neighbours
  if (!already_done)
  {
    for (unsigned int n = 0; n < elem->n_sides(); ++n)
    {
      if (n != entryside)
      {
        const Elem* neigh = elem->neighbor(n);
        if (neigh != NULL)
        {
          int neigh_entry;
          for (unsigned int ns = 0; ns < neigh->n_sides(); ++ns)
          {
            if (neigh->neighbor(ns) == elem)
            {
              neigh_entry = ns;
              break;
            }
          }

          process_element(neigh, neigh_entry, ordered_solutions);
        }

      }
    }
  }
}






void EigenvalueProblem::calculate_dos(void)
{
  if (!get_options().has_submodel("DOS")) return;

  Messages::info("Compute DOS ...");

  ModelOptions& opts = get_options().submodels_begin("DOS")->second;
  
  // create a serial communicator from Device communicator
  libMesh::Parallel::Communicator& comm = get_communicator();
  libMesh::Parallel::Communicator serial_comm;
  comm.split(0,0,serial_comm); 

  delete _energy_mesh;
  _energy_mesh = new Mesh(serial_comm,1);

  double sigma = opts.get_option("gaussian_width", 0.01);

  double emin = opts.get_option("Emin", 0.0);
  double emax = opts.get_option("Emax", 5.0);
  unsigned int num_elem = static_cast<unsigned int>((emax - emin) /
      opts.get_option("dE", 0.001));

  MeshTools::Generation::build_cube (*_energy_mesh,
                                     num_elem, 0, 0,
                                     emin, emax,
                                     0, 0,
                                     0, 0,
                                     EDGE2);


  //
  // The simple approach integrates in k-space with a gaussian weight
  //
  Kspace* kspace_old = _kspace;

  ModelOptions kopts;
  //if (get_options().has_submodel("k-space"))
  //  kopts += get_options().submodels_begin("k-space")->second;

  if (opts.has_submodel("k-space"))
    kopts += opts.submodels_begin("k-space")->second;

  kopts.set_option("mesh_order", "second");
  kopts += parse_kspace_options(kopts);

  _kspace = new Kspace(kopts, serial_comm);


  const MeshBase* kmesh = _kspace->get_k_mesh();
  unsigned int number_of_k_points = kmesh->n_nodes();

  //typedef boost::shared_ptr<eigen_problem_solution> ptr_type;

  //vector<vector<ptr_type>> solutions(number_of_k_points);
  vector<vector<eigen_problem_solution>> solutions(number_of_k_points);

  for (unsigned int i = 0; i < number_of_k_points; i++)
  {
    const Point&  k_point = kmesh->point(i);

    solve_for_kpoint(k_point);
    int number_of_eigs = get_num_states();
    solutions[i].resize(number_of_eigs);

    for (unsigned int j = 0 ; j < number_of_eigs; j++)
      solutions[i][j] = _solution[j];
//      solutions[i][j] = ptr_type(new eigen_problem_solution(_solution[j]));

  }

  // order the solutions according to the bands
  //process_element(kmesh->elem(0), 0, solutions);

  const unsigned int n_energy = _energy_mesh->n_nodes();
  vector<double> dos(n_energy, 0.0);


  double a = 1.0 / (sigma * sqrt(2*M_PI));
  double factor = 1.0 / (2 * M_PI);
  switch (_kspace->dimension())
  {
    case 3:
      a *= factor;
    case 2:
      a *= factor;
    default:
      a *= factor;
  }



  UniquePtr<FEBase> fe(FEBase::build(_kspace->dimension(),
      FEType(SECOND, LAGRANGE)));

  QGauss qrule(_kspace->dimension(),
      static_cast<libMeshEnums::Order>(opts.get_option("integration_order", 17)));
  fe->attach_quadrature_rule(&qrule);

  const std::vector<Real>& JxW = fe->get_JxW();
  //const std::vector<Point>& q_point = fe->get_xyz();
  const std::vector<std::vector<Real> >& phi = fe->get_phi();

  MeshBase::const_element_iterator kit(kmesh->elements_begin());
  const MeshBase::const_element_iterator kend(kmesh->elements_end());
  for ( ; kit != kend; ++kit)
  {
    const Elem* kelem = *kit;

    fe->reinit(kelem);
    unsigned int n_eigs = solutions[kelem->node(0)].size();

    for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
    {
      double w = JxW[qp];

      vector<double> energy(n_eigs, 0.0);
      for (unsigned int i = 0; i < kelem->n_nodes(); ++i)
      {
        unsigned int node = kelem->node(i);
        for (unsigned int n = 0; n < n_eigs; ++n)
        {
          energy[n] += solutions[kelem->node(i)][n].eigen_energy * phi[i][qp];
        }
      }

      for (unsigned int n = 0; n < n_energy; n++)
      {
        double erg =  _energy_mesh->point(n)(0);

        double sum = 0.0;

        for (unsigned int k = 0; k < n_eigs; ++k)
        {
          double ediff = (erg - energy[k]) / sigma;
          double arg = 0.5 * ediff * ediff;
          sum += exp(-arg);
        }

        dos[n] += w * a * sum;
      }
    }

  }

  std::string filename(get_name() + "_dos");

  DataOutput data_output(*_energy_mesh, "dat");
  data_output.set_output_directory(get_output_directory());
  //data_output.set_filename(filename);

  vector<string> names(1, "DOS");
  data_output.write_nodal_data(filename, dos, names);


  // this plots the dispersion
  {
    std::vector<double> results;
    std::vector<std::string> names;

    unsigned int number_of_eigs = solutions[0].size();
    names.resize(number_of_eigs);

    unsigned int number_of_k_points = kmesh->n_nodes();
    results.resize( number_of_eigs * number_of_k_points );


    for (unsigned int i = 0; i < number_of_eigs ; i++)
    {
      std::ostringstream i_str;
      //The states are numbered starting from 0
      i_str << "state_number_" << i;
      names[i] = i_str.str();

      for (unsigned int j = 0; j < number_of_k_points ; j++)
        results[number_of_eigs * j + i] = solutions[j][i].eigen_energy;
    }


    std::string filename(get_name() + "_dispersion");

    DataOutput data_output(*kmesh, "vtk");
    data_output.set_output_directory(get_output_directory());
    //data_output.set_filename(filename);

    data_output.write_nodal_data(filename, results, names);
  }

  _kspace = kspace_old;
}


ID
EigenvalueProblem::do_remember_current_solution(ID id)
{
  map<ID, vector<eigen_problem_solution>>::iterator end(_remembered_sol.end());
  map<ID, vector<eigen_problem_solution>>::iterator it(_remembered_sol.find(id));

  if (it != end)
    it->second = _solution;
  else
  {
    if (_remembered_sol.begin() == end)
      id = 1;
    else
      id = (--end)->first + 1;

    _remembered_sol[id] = _solution;
  }


  return id;
}


void
EigenvalueProblem::do_set_to_remembered_solution(ID id)
{
  map<ID, vector<eigen_problem_solution>>::iterator end(_remembered_sol.end());
  map<ID, vector<eigen_problem_solution>>::iterator it(_remembered_sol.find(id));

  if (it != end)
    _solution = it->second;
}


void
EigenvalueProblem::do_delete_remembered_solution(ID id)
{
  map<ID, vector<eigen_problem_solution>>::iterator end(_remembered_sol.end());
  map<ID, vector<eigen_problem_solution>>::iterator it(_remembered_sol.find(id));

  if (it != end)
    _remembered_sol.erase(it);
}


void EigenvalueProblem::do_plot(void)
{

  SimulationInterface::do_plot();

  if(get_options().has_submodel("Dispersion"))
  {
    Messages m;
    m.info("Compute Dispersion ...");
    m.indent();
  
    ModelOptions::submodel_iterator it(get_options().submodels_begin("Dispersion"));
    const ModelOptions& opts = it->second;
    
    // Back up model kspace
    Kspace* original_kspace = _kspace;

    compute_dispersion(opts);
    
    std::string filename(get_name() + "_dispersion");
    plot_dispersion(filename);

    delete _kspace;
    _kspace = original_kspace;
  }  
  
  if(get_options().has_submodel("BulkDispersion"))
  {
    ModelOptions::submodel_iterator it(get_options().submodels_begin("BulkDispersion"));
    ModelOptions::submodel_iterator end(get_options().submodels_end("BulkDispersion"));
    
    // Back up model kspace
    Kspace* original_kspace = _kspace;
    JobKind original_job = _job;
  

    for( ; it != end ; it++)
    {  
       ModelOptions& opts = it->second;  

       _job = BULKEIGENSTATES;
       vector<double> point(3,0.0);
       opts.get_option("point",point);
       _bulk_point(0) = point[0]; 
       _bulk_point(1) = point[1]; 
       _bulk_point(2) = point[2]; 
       
       ostringstream os; 
       os<<"("<<point[0]<<","<<point[1]<<","<<point[2]<<")";
       
       Messages m;
       m.info("Compute Bulk Dispersion at point "+ os.str() +" ...");
       m.indent();
       opts.set_option("k_space_dimension",3);
   
       Messages::info("Compute dispersion");
       compute_dispersion(opts);
    
       std::string filename(get_name() + "_dispersion_" + os.str() );
       Messages::info("Plot dispersion");
       plot_dispersion(filename);

       delete _kspace;
    }

    _kspace = original_kspace;
    _job = original_job;
  }

  calculate_dos();
}






void
EigenvalueProblem::integrate_density(DofField& density)
{
  // maybe this stuff should be taken from the intial/final state models?
  unsigned int mesh_dim = get_mesh().mesh_dimension();

  bool x_per = get_option("x-periodicity", false);
  bool y_per = get_option("y-periodicity", (mesh_dim < 2));
  bool z_per = get_option("z-periodicity", (mesh_dim < 3));

  ModelOptions kopts;

  if (get_options().has_submodel("k_integration"))
  {
    ModelOptions::submodel_iterator it(get_options().submodels_begin("k_integration"));
    kopts = it->second;
  }
  else
    kopts.set_option("gamma_point_calculation", true);

  kopts.set_option("mesh_units", get_mesh_units());


  int verbose = kopts.get_option("verbose", SimulationOptions::verbose());
  kopts.set_option("verbose", verbose);

  // the normalization volume (in nm)
  double normalization_volume = 1;

  // these are the real space lattice vectors, in nm
  // why pi? Because then the default max k becomes 1 ( = 2*pi/(2*a) ), and
  // k_max can be interpreted in nm^-1
  RealVectorValue a(M_PI, 0.0, 0.0), b(0.0, M_PI, 0.0), c(0.0, 0.0, M_PI);
  auto bbox = get_environment().get_bounding_box();
  if (x_per && (mesh_dim > 0))
  {
    a(0) = (bbox.second(0) - bbox.first(0)) * get_mesh_units() * 1e9;
    normalization_volume *= a(0) * 1e-9;
  }
  if (y_per && (mesh_dim > 1))
  {
    b(1) = (bbox.second(1) - bbox.first(1)) * get_mesh_units() * 1e9;
    normalization_volume *= b(1) * 1e-9;
  }
  if (z_per && (mesh_dim > 2))
  {
    c(2) = (bbox.second(2) - bbox.first(2)) * get_mesh_units() * 1e9;
    normalization_volume *= c(2) * 1e-9;
  }

  // if there is an atomistic structure, we can take the lattice vectors from it
  // (they come in Angstrom!)
  if (get_atomistic_structure() != NULL)
  {
    const AtomisticStructure* as = get_atomistic_structure();
    as->get_lattice_vectors(a, b, c);
    // scale to nm
    a *= 0.1;
    b *= 0.1;
    c *= 0.1;

    normalization_volume = 1;
    if ((x_per = as->is_periodic(0)))
      normalization_volume *= a(0);
    if ((y_per = as->is_periodic(1)))
      normalization_volume *= b(1);
    if ((z_per = as->is_periodic(2)))
      normalization_volume *= c(2);
  }

  kopts.set_option("normalization_volume", normalization_volume);


  unsigned int k_dim = 0;
  if (x_per)
    k_dim++;
  if (y_per)
    k_dim++;
  if (z_per)
    k_dim++;

  //k_dim = min(k_dim, 3u);

  if (kopts.find_option("k_space_dimension"))
    k_dim = kopts.get_option("k_space_dimension", k_dim);

  kopts.set_option("k_space_dimension", k_dim);

  switch (k_dim)
  {
    case 1:
      if (mesh_dim == 3)
      {
        if (x_per)
          c = a;
        else if (y_per)
          c = b;
      }
      kopts.set_option("r1", c);
      break;

    case 2:
      if (mesh_dim == 3)
      {
        if (!z_per)
        {
          c = b;
          b = a;
        }
        else if (x_per)
        {
          b = a;
        }
      }
      else if (mesh_dim == 2)
      {
        if (x_per)
          b = a;
      }
      kopts.set_option("r1", b);
      kopts.set_option("r2", c);
      break;

    case 3:
      kopts.set_option("r1", a);
      kopts.set_option("r2", b);
      kopts.set_option("r3", c);
      break;

    default:
      break;
  }


  Messages m;
  m.info("Setting up k-space integration");
  m.indent();

  KspaceIntegration* kint = KspaceIntegration::create(this,
      &EigenvalueProblem::calculate_density_at_k, kopts, get_communicator(), get_solver_communicator());

  if (kint == NULL)
    throw InitFailedException("Could not create k-integration for density calculation");

  kint->init();

  kint->solve();

  density = kint->get_solution();
}


void
EigenvalueProblem::calculate_density_at_k(const Point& k_point,
        DofField& density, double& )
{
  solve_for_kpoint(k_point);
  do_calculate_density_at_k(density);
}


 
void EigenvalueProblem::solve_for_kpoint(const Point& kpoint)
{
  const Point oldk(_k_vector[0], _k_vector[1], _k_vector[2]);
  set_k_point(kpoint);
  do_solve_for_kpoint(kpoint);
  ostringstream os;
  if (get_option("plot_at_every_k", false))
  {
    os << "k(" << _k_vector[0] << "," <<  _k_vector[1] << "," <<
        _k_vector[2] << ")";
    TiberCad::prepend_to_filename_suffix(os.str());
    this->plot_meshdata();
    this->plot_globaldata();
    this->plot_atomisticdata();
    TiberCad::drop_first_filename_suffix();
  }
  set_k_point(oldk);
  k_is_old();
}


void
EigenvalueProblem::do_solve_for_kpoint(const Point& )
{
  solve();
}


void EigenvalueProblem::get_eigenvalues(const std::string& particle, 
					std::vector<double>& values) const
{

  unsigned int n = _solution.size();
  unsigned int num_st = 0;
  values.reserve(n);

  for (unsigned int i = 0; i < n; i++)
  {
    if(particle.empty() || (_solution[i].particle == particle))
    {  
      num_st++;

      values.push_back( _solution[i].eigen_energy ); 
    }
  }  
 
  values.resize(num_st);

}

unsigned int EigenvalueProblem::get_num_states(void) const
{
  return _solution.size();
}


unsigned int EigenvalueProblem::get_num_states(const std::string& particle) const
{
  unsigned int num_i_states = 0;
  for(unsigned int i=0; i<_solution.size(); i++)
  {
    if(_solution[i].particle == particle) num_i_states++;  
  }
  
  return num_i_states;
}

std::vector<ID>
EigenvalueProblem::get_state_indices(const std::string& particle) const
{
  unsigned int num = get_num_states(particle);	
  std::vector<ID> result(num, 0);

  unsigned int num_st=0;
  for(unsigned int i=0; i<_solution.size(); i++)
  {
    if(_solution[i].particle == particle)
      { result[num_st]=i; num_st++; }
  }
  
  return result;
}


void EigenvalueProblem::get_populations(const std::string& particle, 
					std::vector<double>& values) const
{
 
  unsigned int n = _solution.size();
  unsigned int num_st = 0;
  values.reserve(n);
 
  for (unsigned int i = 0; i < n; i++)
  {
    if(_solution[i].particle == particle)
    {  
      num_st++;

      if(_solution[i].statistics == "Fermi")
      {      
	double val = Fermi(_solution[i].eigen_energy, _solution[i].electro_chem_pot, 
		       _solution[i].temperature);

	if(particle == "el" || particle == "electron")
	{
	  values.push_back(val);	  
	}	
	
	if(particle == "hl" || particle == "hole")
	{
	  values.push_back(1-val);	  
	}

      }
      else
      {
	double val = Bose(_solution[i].eigen_energy, _solution[i].electro_chem_pot, 
		       _solution[i].temperature);

	values.push_back(val);	
	
      }

    }
     
  }

  values.resize(num_st);
 
} 
 
double  EigenvalueProblem::get_population(int i) const
{

  double val = 0.0;

  if(_solution[i].statistics == "Fermi")
  {        
    val = Fermi(_solution[i].eigen_energy, _solution[i].electro_chem_pot,
        _solution[i].temperature);


    if(_solution[i].particle == "hl" || _solution[i].particle == "hole")
    {
      val = 1 - val;
    }

  }
  else
  {
    val = Bose(_solution[i].eigen_energy, _solution[i].electro_chem_pot,
        _solution[i].temperature);

  }

  return val;
}




double  EigenvalueProblem::Fermi(double Energy, double Fermi_energy, double Temperature) const
{
  double T_EV = Temperature * Constants::k_Boltzmann;
  double exp_arg =  (Energy - Fermi_energy)/T_EV;
  
  double occupation;
  
  if (exp_arg > 35) 
    occupation = std::exp(-exp_arg);
  else
    occupation = 1.0/(std::exp(exp_arg) + 1.0);
  
  return occupation;

}

double  EigenvalueProblem::Bose(double Energy, double electro_chem_pot, double Temperature) const
{
  double T_EV = Temperature * Constants::k_Boltzmann;
  double exp_arg =  (Energy - electro_chem_pot)/T_EV;
  
  double bose;
  
  if (exp_arg > 35) 
    bose = std::exp(-exp_arg);
  else
    bose = 1.0/(std::exp(exp_arg) - 1.0);
  
  return bose;

}


void EigenvalueProblem::write_states(void) const
{

  int num_st=_solution.size();

  Messages::newline();
  Messages::info("#  type   level    stat.     pot.       pop.");


  for(int i=0; i< num_st; i++)
  {
    ostringstream os;
    os << i << " " << _solution[i].particle << " " << std::setprecision(6)
	     << _solution[i].eigen_energy << " " <<_solution[i].statistics
	     << " " <<std::setw(10) << _solution[i].electro_chem_pot
	     << " " <<std::setw(10) << get_population(i);
    Messages::info(os.str());
  }
  Messages::newline();

}



  
void EigenvalueProblem::copy_H_to_solver( )
{
  do_copy_H_to_solver();
}

void EigenvalueProblem::copy_S_to_solver( )
{
  do_copy_S_to_solver();
}





Complex
EigenvalueProblem::scalar_product(const eigen_problem_solution& a,
                                  const eigen_problem_solution& b) const
{
  return scalar_product(a.eigen_vector, b.eigen_vector);
}



Complex
EigenvalueProblem::scalar_product(const vector<Complex>& a,
                                  const vector<Complex>& b) const
{
  Complex sprod(0,0);

  if (a.size() == b.size())
  {
    size_t length = a.size();
    for (size_t i = 0; i < length; ++i)
      sprod += conj(a[i]) * b[i];
  }

  return sprod;
}
