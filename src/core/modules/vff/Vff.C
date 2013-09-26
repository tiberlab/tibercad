// $Id$

#include "Vff.h"
#include "VffModel.h"
#include "Messages.h"
#include "AtomisticStructure.h"
#include "Specie.h"
#include "OptGpl.h"
#include "Atom.h"
#include "RuntimeException.h"
#include "TiberModule.h"
#include "point.h"
#include <cmath>
#include "mesh.h"


Vff*
Vff::_this = NULL;


Vff::Options::Options(void)
: boundary_conditions("free_standing"),
  substrate_plane("z"),
  substrate_tol(1.0),
  absolute_tolerance(1e-3),
  method("cg")
{
}

Vff::Vff(const ModelOptions& options) :
              SimulationInterface(options),
              _has_strain_tensor(false)
{
  // there's nothing to be done
}

Vff::~Vff(void)
{
  // there's nothing to be done
}

Vff*
Vff::create(const ModelOptions& options)
{
  return new Vff(options);
}

PhysicalModel*
Vff::create_bulk_model(const ModelOptions& options,
    const Material* mat) const
{

  return VffModel::create(mat, options);
}

void
Vff::parse_options()
{
  const ModelOptions& opts = get_options();
  Options& myopts = get_my_options();

   myopts.boundary_conditions = opts.get_option("boundary_conditions", "free_standing");
   myopts.substrate_plane = opts.get_option("substrate_plane", "z");
   myopts.substrate_tol = opts.get_option("substrate_tol", 1.0);
   myopts.substrate_updown = opts.get_option("substrate_updown", false);

   //Solver options
   myopts.method = get_solver_options().get_option("method", "cg");
   myopts.absolute_tolerance = get_solver_options().get_option("absolute_tolerance", 1e-3);
}

void
Vff::do_init()
{
  Messages::info("Initializing VFF module");

  if (get_atomistic_structure()==NULL)
    throw InitFailedException("VFF: could not find atomistic structure");

  Messages::debug("Setting boundary conditions");
  check_structure();
  parse_options();
  declare_solution(StrainNodes, NTUPLE, NODES,"unitless",6);
  declare_solution(StrainCells, NTUPLE, CELL,"unitless",6);
}

void
Vff::do_solve(void)
{
  Messages::debug("Initializing VFF Module");
   //Set_boundary temporary here otherwise the dof are fixed before Macrostrain runs
   set_boundary();

   //Build preliminar information
   resize_parameters();
   build_parameters();
   set_coords();

   Messages::debug("Starting VFF Conjugate Gradient");
   //Calling optimization
   optimize();

   //Apply displacement to AtomisticStructure instance
   displace_atoms();

}

void
Vff::set_coords(void)
{
  _n_atoms = get_atomistic_structure()->get_N_without_H();
  unsigned int n_all_atoms = get_atomistic_structure()->get_N_atoms();
  unsigned int n_atoms = _n_atoms;
  std::vector<double>& coords = get_coords();

  coords.resize(n_atoms * 3, 0.0);
  _initial_coords.resize(n_all_atoms * 3, 0.0);
  //Put all coords in a temporary 1D array
  for (unsigned int i = 0; i < n_atoms; i ++)
    {
      unsigned int i_start = i * 3;
      coords[i_start] =
          get_atomistic_structure()->get_structure_atoms()[i].get_position(0);
      coords[i_start + 1] =
          get_atomistic_structure()->get_structure_atoms()[i].get_position(1);
      coords[i_start + 2] =
          get_atomistic_structure()->get_structure_atoms()[i].get_position(2);
    }

  for (unsigned int i = 0; i < n_all_atoms; i ++)
    {
      unsigned int i_start = i * 3;
      _initial_coords[i_start] =
          get_atomistic_structure()->get_structure_atoms()[i].get_position(0);
      _initial_coords[i_start + 1] =
          get_atomistic_structure()->get_structure_atoms()[i].get_position(1);
      _initial_coords[i_start + 2] =
          get_atomistic_structure()->get_structure_atoms()[i].get_position(2);
    }

}


void
Vff::displace_atoms(void)
{
  std::vector<Atom>& atoms = get_atomistic_structure()->get_structure_atoms();
  std::vector<unsigned int>& free_atoms = get_free_atoms();

  //Displace hydrogens according to bonded atom
  const Bondmap& bondmap = get_atomistic_structure()->get_bond_map();
  double x_t, y_t, z_t;
  unsigned int n_atoms = get_atomistic_structure()->get_N_without_H();
  unsigned int n_all_atoms = get_atomistic_structure()->get_N_atoms();

  //Displace other atoms
  for (unsigned int i = 0; i < _n_free_atoms; i++)
    {
      unsigned int j = free_atoms[i];
      atoms[j].set_position(0, _dof[i * 3 + 0]);
      atoms[j].set_position(1, _dof[i * 3 + 1]);
      atoms[j].set_position(2, _dof[i * 3 + 2]);
    }

  for (unsigned int i = n_atoms + 1; i < n_all_atoms; i++)
      {
      unsigned int n_bonds = bondmap[i].size();
      if (n_bonds != 1)
        throw RuntimeException("Structure is badly defined. Passivation atom has multiple bonds");
      unsigned int j = bondmap[i][0];
      //Calculate translation
      x_t = _initial_coords[j * 3 + 0] - atoms[j].get_position(0);
      y_t = _initial_coords[j * 3 + 1] - atoms[j].get_position(1);
      z_t = _initial_coords[j * 3 + 2] - atoms[j].get_position(2);
      //Translate hydrogen accordingly
      atoms[i].set_position(0, _initial_coords[i * 3 + 0] + x_t);
      atoms[i].set_position(1, _initial_coords[i * 3 + 1] + y_t);
      atoms[i].set_position(2, _initial_coords[i * 3 + 2] + z_t);
      }

  _initial_coords.resize(0);
  get_coords().resize(0);
  _dof.resize(0);

  get_atomistic_structure()->print_structure("strained_VFF.xyz");
}

void
Vff::set_boundary(void)
{
  //WARNING: passivation hydrogens are always left out!!!

  int n_atoms = get_atomistic_structure()->get_N_without_H();
  const Bondmap& bondmap = get_atomistic_structure()->get_bond_map();

  std::vector<unsigned int>& free_atoms = get_free_atoms();
  std::vector<Atom>& atoms = get_atomistic_structure()->get_structure_atoms();

  free_atoms.reserve(n_atoms);

  //Assign free atoms indexes according to selected boundary conditions
  if (get_my_options().boundary_conditions == "free_standing")
    {
      for (unsigned int i = 0; i < n_atoms; i++) free_atoms.push_back(i);
      _n_free_atoms = n_atoms;
    }

  if  (get_my_options().boundary_conditions == "all_around_1atom")
    {
      _n_free_atoms = 0;
      for (unsigned int i = 0; i < n_atoms; i++)
        {
          //If bonded to 4 non-hydrogen atoms, then it's a free atom
          if ((bondmap[i].size() == 4) &&
              (get_atomistic_structure()->get_specie(bondmap[i][0]) != Specie::H ) &&
              (get_atomistic_structure()->get_specie(bondmap[i][1]) != Specie::H ) &&
              (get_atomistic_structure()->get_specie(bondmap[i][2]) != Specie::H ) &&
              (get_atomistic_structure()->get_specie(bondmap[i][3]) != Specie::H ))
            {
              free_atoms.push_back(i);
              _n_free_atoms += 1;
            }
        }
    }



  if  (get_my_options().boundary_conditions == "all_around")
    {

      _n_free_atoms = 0;
      for (unsigned int i = 0; i < n_atoms; i++)
        {
          bool free = true;
          for (unsigned int counter_j = 0; counter_j < bondmap[i].size(); counter_j++)
            {
              unsigned int j = bondmap[i][counter_j];
              for (unsigned int counter_k = 0; counter_k < bondmap[j].size(); counter_k++)
                {
                  unsigned int k = bondmap[j][counter_k];
                  if ((bondmap[k].size() != 4) ||
                      (get_atomistic_structure()->get_specie(bondmap[k][0]) == Specie::H ) ||
                      (get_atomistic_structure()->get_specie(bondmap[k][1]) == Specie::H ) ||
                      (get_atomistic_structure()->get_specie(bondmap[k][2]) == Specie::H ) ||
                      (get_atomistic_structure()->get_specie(bondmap[k][3]) == Specie::H ))
                    {
                      free = false;
                    }
                }
            }

          //If bonded to 4 non-hydrogen atoms, then it's a free atom
          if (free)
            {
              free_atoms.push_back(i);
              _n_free_atoms += 1;
            }
        }
    }
  if (get_my_options().boundary_conditions == "substrate")
   {
     double tol = get_my_options().substrate_tol;
     unsigned int plane = 0;
     if (get_my_options().substrate_plane == "z")
         plane = 2;
     if (get_my_options().substrate_plane == "y")
             plane = 1;
     if (get_my_options().substrate_plane == "x")
             plane = 0;
     double min_coord =
         get_atomistic_structure()->get_structure_atoms()[0].get_position(plane);
     double max_coord =
         get_atomistic_structure()->get_structure_atoms()[0].get_position(plane);
     bool condition = false;
     for (unsigned int i = 0; i < n_atoms; i++)
     {
       double coord =
           get_atomistic_structure()->get_structure_atoms()[i].get_position(plane);
       if (coord < min_coord)
           min_coord = coord;
       if (coord > max_coord)
         max_coord = coord;
     }
     for (unsigned int i = 0; i < n_atoms; i++)
          {
            double coord =
                get_atomistic_structure()->get_structure_atoms()[i].get_position(plane);
            if (!get_my_options().substrate_updown)
              condition = (coord > min_coord + tol);
            if (get_my_options().substrate_updown)
              condition = (coord > min_coord + tol) && (coord < max_coord - tol);
            if (condition)
            {
              free_atoms.push_back(i);
            }
          }
   }

  _n_free_atoms = free_atoms.size();

  //Common operations
  _n_dof = _n_free_atoms * 3;
  _dof.resize(_n_dof, 0.0);
  for (unsigned int i = 0; i < _n_free_atoms; i++)
    {
      unsigned int j = free_atoms[i];
      _dof[i * 3 + 0] = atoms[j].get_position(0);
      _dof[i * 3 + 1] = atoms[j].get_position(1);
      _dof[i * 3 + 2] = atoms[j].get_position(2);
    }

  //Write some informations
  Messages::info("VFF boundary conditions set up");
  std::cout << "Number of free atoms: " << _n_free_atoms << std::endl;

}

void
Vff::check_structure(void)
{

  Specie::Type h_specie = Specie::H;

  //Check that no hydrogen is here
  for (unsigned int i = 0; i < get_atomistic_structure()->get_N_without_H(); i++)
    {
      if (get_atomistic_structure()->get_specie(i) == Specie::H)
        throw InitFailedException("VFF: hydrogens with wrong array index the structure");
    }


}

void
Vff::resize_parameters(void)
{
  unsigned int n_max_neighbors = 4;
  int n_atoms = get_atomistic_structure()->get_N_without_H();

  _alpha.resize(n_atoms);
  _d.resize(n_atoms);
  for (unsigned int i = 0; i < n_atoms; i++)
    {
      _alpha[i].resize(n_max_neighbors, 0.0);
      _d[i].resize(n_max_neighbors, 0.0);
    }

  _beta.resize(n_atoms);
  _teta.resize(n_atoms);
  for (unsigned int i = 0; i < n_atoms; i++)
    {
      _beta[i].resize(n_max_neighbors);
      _teta[i].resize(n_max_neighbors);
      for (unsigned j = 0; j < n_max_neighbors; j++)
        {
          _beta[i][j].resize(n_max_neighbors, 0.0);
          _teta[i][j].resize(n_max_neighbors, 0.0);
        }
    }
}

void
Vff::build_parameters(void)
{
  int n_atoms = get_atomistic_structure()->get_N_without_H();
  const Bondmap& bondmap = get_atomistic_structure()->get_bond_map();
  VffModel* pm_a = NULL;
  VffModel* pm_b = NULL;

  bool parent(get_atomistic_structure()->is_random_alloy());

  for (unsigned int i = 0; i < n_atoms; i++)
    {
      const Atom& atm_i = get_atomistic_structure()->get_structure_atom(i);
      unsigned int n_bonds = bondmap[i].size();
      for (unsigned int counter_j = 0; counter_j < n_bonds; counter_j++)
        {
          unsigned int j = bondmap[i][counter_j];
          //You don't want informations for passivation bonds
          if (get_atomistic_structure()->get_specie(j) != Specie::H)
            {
              const Atom& atm_j = get_atomistic_structure()->get_structure_atom(j);
              pm_a = get_bulk_model<VffModel>(atm_i, atm_j, parent);
              _alpha[i][counter_j] = pm_a->get_alpha(atm_i, atm_j);
              _d[i][counter_j] = pm_a->get_d(atm_i, atm_j);

              for (unsigned int counter_k = 0; counter_k < n_bonds; counter_k++)
                {
                  unsigned int k = bondmap[i][counter_k];
                  const Atom& atm_k = get_atomistic_structure()->get_structure_atom(k);
                  if (get_atomistic_structure()->get_specie(k) != Specie::H)
                    {
                      pm_a = get_bulk_model<VffModel>(atm_i, atm_j, parent);
                      pm_b = get_bulk_model<VffModel>(atm_i, atm_k, parent);

                      _teta[i][counter_j][counter_k] =
                          (pm_a->get_costeta(atm_i, atm_j, atm_k) +
                              pm_b->get_costeta(atm_i, atm_j, atm_k)) / 2.0;
                      _beta[i][counter_j][counter_k] = sqrt(
                          pm_a->get_beta(atm_i, atm_j, atm_k) *
                          pm_b->get_beta(atm_i, atm_j, atm_k));
                    }

                }

            }
        }
    }
}


double
Vff::keating_potential(void)
{
  double u = 0.0;
  //Atoms to be considered in keating potential (H passivation not included)
  int n_atoms = _n_atoms;
  const Bondmap& bondmap = get_atomistic_structure()->get_bond_map();
  const std::vector<std::vector<Tensor1> > translation =
      get_atomistic_structure()->get_neighbor_translation();
  std::vector<double>& coords = get_coords();


  //Substitute degree of freedom (they can change during optimization) in local coords
  for (unsigned int i = 0; i < _n_free_atoms; i ++)
    {
      unsigned int free_index = i * 3;
      unsigned int index = _free_atoms[i] * 3;
      coords[index] = _dof[free_index]; coords[index + 1] = _dof[free_index + 1];
      coords[index + 2] = _dof[free_index + 2];

    }

  //Calculate the potential
  //NOTE: i, j and k denote an atom index. counter_j and counter_k keep track of the index of neighbor (0, 1, 2, 3)
  for (unsigned int i = 0; i < n_atoms; i++)
    {


      unsigned int n_bonds = bondmap[i].size();

      for (unsigned int counter_j =0; counter_j < n_bonds; counter_j++)
        {

          unsigned int j = bondmap[i][counter_j];
          double t_j_x = translation[i][counter_j](1);
          double t_j_y = translation[i][counter_j](2);
          double t_j_z = translation[i][counter_j](3);
          // t_j_x = 0.0;
          // t_j_y = 0.0;
          // t_j_z = 0.0;
          //Hydrogen bonds must not be included (they don't have reference distance and angles)
          if (get_atomistic_structure()->get_structure_atoms()[j].get_specie() != Specie::H)
            {


              unsigned int i_start = i * 3; unsigned int j_start = j * 3;

              double prefactor = (_alpha[i][counter_j] * 3.0) / (16.0 * _d[i][counter_j] * _d[i][counter_j]);
              double x_ij = coords[i_start] - coords[j_start] - t_j_x;
              double y_ij = coords[i_start + 1] - coords[j_start + 1] - t_j_y;
              double z_ij = coords[i_start + 2] - coords[j_start + 2] - t_j_z;

              double bond_stretching = prefactor *
                  pow((x_ij * x_ij + y_ij * y_ij + z_ij * z_ij - _d[i][counter_j] * _d[i][counter_j]),2);

              u += bond_stretching;

              for (unsigned int counter_k = 0; counter_k < n_bonds; counter_k++)
                {
                  if (counter_k != counter_j)
                    {

                      unsigned int k = bondmap[i][counter_k];
                      double t_k_x = translation[i][counter_k](1);
                      double t_k_y = translation[i][counter_k](2);
                      double t_k_z = translation[i][counter_k](3);
                      // t_k_x = 0.0;
                      // t_k_y = 0.0;
                      // t_k_z = 0.0;
                      //Hydrogen bonds must not be included (they don't have reference distance and angles)
                      if (get_atomistic_structure()->get_structure_atoms()[k].get_specie() != Specie::H)
                        {
                          unsigned int k_start = k * 3;
                          double x_ik = coords[i_start] - coords[k_start] - t_k_x;
                          double y_ik = coords[i_start + 1] - coords[k_start + 1] - t_k_y;
                          double z_ik = coords[i_start + 2] - coords[k_start + 2] - t_k_z;

                          prefactor = (_beta[i][counter_j][counter_k] * 3.0) / (16.0 * _d[i][counter_j] * _d[i][counter_k]);

                          double bond_bending = prefactor *
                              pow((x_ij * x_ik + y_ij * y_ik + z_ij * z_ik - _d[i][counter_j] * _d[i][counter_k] * _teta[i][counter_j][counter_k]),2);

                          u += bond_bending;

                        }
                    }
                }

            }
        }

    }
  std::cout << "potential calculated " << u << std::endl;

  return u;


}

double
Vff::keating_potential(double* x, int n)
{
  if (n != _n_dof)
    throw InitFailedException("Cannot calculate keating potential. Number of coordinates and degree of freedoms mismatch");

  for (unsigned int i = 0; i < n; i++)
    {
      _dof[i] = x[i];
    }

  return keating_potential();

}

void
Vff::keating_gradient(double* grad, double* x, int n)
{
  if (n != _n_dof)
    throw InitFailedException("Cannot calculate keating potential. Number of coordinates and degree of freedoms mismatch");

  for (unsigned int i = 0; i < n; i++)
    {
      _dof[i] = x[i];
    }

  std::vector<double> grad_vec(keating_gradient());

  for (unsigned int i = 0; i < n; i++)
    {
      grad[i] = grad_vec[i];
    }
}

double
Vff::keating_pot_grad(double* grad, double* x, int n)
{
  keating_gradient(grad, x, n);
  return keating_potential(x,n);
}


std::vector<double>
Vff::keating_gradient(void)
{

  //Atoms to be considered in gradient (H passivation not included)
  int n_atoms = _n_atoms;
  const Bondmap& bondmap = get_atomistic_structure()->get_bond_map();
  std::vector<double> coords = get_coords();
  const std::vector<std::vector<Tensor1> > translation =
      get_atomistic_structure()->get_neighbor_translation();

  //Note: grad_all is introduced as it's easier coding the gradient in a general way
  //considering all the coordinates. Then we only take the elements corresponding to degrees
  //of freedom
  std::vector<double> grad_all(n_atoms * 3, 0.0);
  std::vector<double> grad_dof(_n_free_atoms * 3, 0.0);


  //Substitute degree of freedom (they can change during optimization) in local coords
  for (unsigned int i = 0; i < _n_free_atoms; i ++)
    {
      unsigned int free_index = i * 3;
      unsigned int index = _free_atoms[i] * 3;
      coords[index] = _dof[free_index]; coords[index + 1] = _dof[free_index + 1];
      coords[index + 2] = _dof[free_index + 2];
    }


  //Calculate the gradient
  //NOTE: i, j and k denote an atom index. counter_j and counter_k keep track of the index of neighbor (0, 1, 2, 3)
  for (unsigned int i = 0; i < n_atoms; i++)
    {

      unsigned int n_bonds = bondmap[i].size();

      for (unsigned int counter_j =0; counter_j < n_bonds; counter_j++)
        {

          unsigned int j = bondmap[i][counter_j];
          double t_j_x = translation[i][counter_j](1);
          double t_j_y = translation[i][counter_j](2);
          double t_j_z = translation[i][counter_j](3);
//           t_j_x = 0.0;
//           t_j_y = 0.0;
//           t_j_z = 0.0;
          //Hydrogen bonds must not be included (they don't have reference distance and angles)
          if (get_atomistic_structure()->get_structure_atoms()[j].get_specie() != Specie::H)
            {

              unsigned int i_start = i * 3; unsigned int j_start = j * 3;

              double prefactor = (_alpha[i][counter_j] * 3.0) / (16.0 * _d[i][counter_j] * _d[i][counter_j]);
              double x_ij = coords[i_start] - coords[j_start] - t_j_x;
              double y_ij = coords[i_start + 1] - coords[j_start + 1] - t_j_y;
              double z_ij = coords[i_start + 2] - coords[j_start + 2] - t_j_z;

              double common = (x_ij * x_ij + y_ij * y_ij + z_ij * z_ij - _d[i][counter_j] * _d[i][counter_j]) * prefactor * 4.0;

              double grad_x = common * x_ij;
              double grad_y = common * y_ij;
              double grad_z = common * z_ij;

              grad_all[i_start] += grad_x; grad_all[i_start + 1] += grad_y; grad_all[i_start + 2] += grad_z;
              grad_all[j_start] -= grad_x; grad_all[j_start + 1] -= grad_y; grad_all[j_start + 2] -= grad_z;

              for (unsigned int counter_k = 0; counter_k < n_bonds; counter_k++)
                {
                  if (counter_k != counter_j)
                    {
                      unsigned int k = bondmap[i][counter_k];
                      double t_k_x = translation[i][counter_k](1);
                      double t_k_y = translation[i][counter_k](2);
                      double t_k_z = translation[i][counter_k](3);
//                       t_k_x = 0.0;
//                       t_k_y = 0.0;
//                       t_k_z = 0.0;
                      //Hydrogen bonds must not be included (they don't have reference distance and angles)
                      if (get_atomistic_structure()->get_structure_atoms()[k].get_specie() != Specie::H)
                        {
                          unsigned int k_start = k * 3;
                          double x_ik = coords[i_start] - coords[k_start] - t_k_x;
                          double y_ik = coords[i_start + 1] - coords[k_start + 1] - t_k_y;
                          double z_ik = coords[i_start + 2] - coords[k_start + 2] - t_k_z;

                          prefactor = (_beta[i][counter_j][counter_k] * 3.0) / (16.0 * _d[i][counter_j] * _d[i][counter_k]);

                          common = (x_ij * x_ik + y_ij * y_ik + z_ij * z_ik -
                              _d[i][counter_j] * _d[i][counter_k] * _teta[i][counter_j][counter_k]) * 2.0 * prefactor;

                          grad_all[i_start] += common * (x_ik + x_ij); grad_all[i_start + 1] += common * (y_ik + y_ij);
                          grad_all[i_start + 2] += common * (z_ik + z_ij);

                          grad_all[j_start] -= common * (x_ik); grad_all[j_start + 1] -= common * (y_ik);
                          grad_all[j_start + 2] -= common * (z_ik);

                          grad_all[k_start] -= common * (x_ij); grad_all[k_start + 1] -= common * (y_ij);
                          grad_all[k_start + 2] -= common * (z_ij);

                        }
                    }
                }

            }
        }

    }

  //Insert grad_dof values in right place
  for (unsigned int i = 0; i < _n_free_atoms; i ++)
    {
      unsigned int free_index = i * 3;
      unsigned int index = _free_atoms[i] * 3;
      grad_dof[free_index] = grad_all[index]; grad_dof[free_index + 1] = grad_all[index + 1];
      grad_dof[free_index + 2] = grad_all[index + 2];

    }


  return grad_dof;

}


void
Vff::optimize(void)
{
  if (get_my_options().method == "cg") 
  {
    //Convert absolute tolerance from eV/A to internal units 10^-20J/A
    double tol = get_my_options().absolute_tolerance * 16.0; 
    OptGpl solver(*this);
    solver.solve(tol);
  }
  else
    Messages::error("Error in VFF: solver method does not exist.");

}


void
Vff::get_solution_secure(const Elem* elem,
    std::map<ID, std::vector<double> >& values,
    const std::vector<Point>& p)
{
//TODO: these methods are experimental, their only purpose is to test different
//techniques, they should not be trusted without speaking with the developer
//
//
if (!_has_strain_tensor)
{
   std::cout << "Invoking StrainLattice " << std::endl;
   _strain.init(get_atomistic_structure());
   std::cout << "Calculating strain tensor "  << std::endl;
   _strain.do_solve();
   std::cout << "Strain tensor available " << std::endl;
   _has_strain_tensor = true;
}

unsigned int np = p.size();
double scale = get_atomistic_structure()->get_scale();
unsigned int dim = get_mesh().mesh_dimension();
Point phys_p;
if (values.count(StrainNodes))
{
for (unsigned int n = 0; n < np; n++)
  {
  if (dim == 1)
   phys_p = FE< 1, libMeshEnums::LAGRANGE>::map(elem, p[n]);
  if (dim == 2)
   phys_p = FE< 2, libMeshEnums::LAGRANGE>::map(elem, p[n]);
  if (dim == 3)
   phys_p = FE< 3, libMeshEnums::LAGRANGE>::map(elem, p[n]);
   double min = 1e5;
   Tensor2Gen sol(0);
   for (unsigned int i=0; i < _strain.get_solution().size(); i++)
   {
     double distance =
     (_strain.get_solution()[i].atom_p->get_position()(0) - phys_p(0)*scale) *
     (_strain.get_solution()[i].atom_p->get_position()(0) - phys_p(0)*scale) +
     (_strain.get_solution()[i].atom_p->get_position()(1) - phys_p(1)*scale) *
     (_strain.get_solution()[i].atom_p->get_position()(1) - phys_p(1)*scale) +
     (_strain.get_solution()[i].atom_p->get_position()(2) - phys_p(2)*scale) *
     (_strain.get_solution()[i].atom_p->get_position()(2) - phys_p(2)*scale);
     if (distance<min) 
     {sol = _strain.get_solution()[i].tensor;min=distance;}     
   }
   values[StrainNodes][6*n]=sol(1,1);
   values[StrainNodes][6*n+1]=sol(2,2);
   values[StrainNodes][6*n+2]=sol(3,3);
   values[StrainNodes][6*n+3]=sol(1,2);
   values[StrainNodes][6*n+4]=sol(1,3);
   values[StrainNodes][6*n+5]=sol(2,3);
}
}

//Element solution, take the average of all the tensors included in the 
//element, if not any take the nearest
if (values.count(StrainCells))
{

  unsigned int contribs = 0;
  
  Tensor2Gen sol(0);
   phys_p = elem->centroid();
   double min = 1e5;

   for (unsigned int i=0; i < _strain.get_solution().size(); i++)
   { 
    if (_strain.get_solution()[i].atom_p->get_elem() == elem)
    {
      contribs +=1;
      sol += _strain.get_solution()[i].tensor;
      //std::cout << "from element " << sol << " contribs " << contribs << std::endl;
    }
   }
   if (contribs == 0)
   {
     contribs = 1;
   for (unsigned int i=0; i < _strain.get_solution().size(); i++)
   { 
     double distance =
     (_strain.get_solution()[i].atom_p->get_position()(0) - phys_p(0)*scale) *
     (_strain.get_solution()[i].atom_p->get_position()(0) - phys_p(0)*scale) +
     (_strain.get_solution()[i].atom_p->get_position()(1) - phys_p(1)*scale) *
     (_strain.get_solution()[i].atom_p->get_position()(1) - phys_p(1)*scale) +
     (_strain.get_solution()[i].atom_p->get_position()(2) - phys_p(2)*scale) *
     (_strain.get_solution()[i].atom_p->get_position()(2) - phys_p(2)*scale);
     if (distance<min) 
     {sol = _strain.get_solution()[i].tensor;min=distance;}
  }
   }
   sol = sol / contribs;
   //std::cout << "sol " <<  sol << " contribs " << contribs <<std::endl;
for (unsigned int n = 0; n < np; n++)
{
   values[StrainCells][0]=sol(1,1);
   values[StrainCells][1]=sol(2,2);
   values[StrainCells][2]=sol(3,3);
   values[StrainCells][3]=sol(1,2);
   values[StrainCells][4]=sol(1,3);
   values[StrainCells][5]=sol(2,3);
   //std::cout << " values " <<  values[Strain2][0] << std::endl;
}
}

}
