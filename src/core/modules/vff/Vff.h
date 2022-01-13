
#ifndef _VFF_H_
#define _VFF_H_

#include "SimulationInterface.h"
#include "StrainLattice.h"



class TBDLLOCAL Vff : public SimulationInterface
{

public:

  //! These are the known solution variables
  /*!
   * This is an enum, but we use the string representation of
   * the enum values to refer to solutions for plotting or
   * for data exchange with other modules.
   *
   * \note Do \em not use (\c INVALID_ID - 1) or the strings \c RegionIDs
   * or \c materials as they are used to plot the materials/region IDs.
   *
   * \note The name "all" is used to plot all solutions
   */
  enum Solutions
  {
    StrainNodes,
    StrainCells,
    Displacement
  };

  //! Destructor
  /*!
   * We do not declare it virtual here, as we will not allow
   * to derive from this class anyway.
   */
  ~Vff(void);

  //! We need a public static creator function
  static Vff* create(const ModelOptions& options);

  void keating_gradient(double* grad, double* x, int n);

  double keating_pot_grad(double* grad, double* x, int n);

  //! Calculate Keating potential
  /*!
   * Keating potential is calculated. Resulting potential is in 10^-20 J (N/m * A^2).
   * Potential in eV can be obtained multiplying by 16
   */
  double keating_potential(double* x, int n);

  int get_n_dof(void) const;

  const std::vector<double>& get_dof(void) const;

protected:

  //! The initialization
  virtual void do_init(void);


  //! Parse the options from the input file
  virtual void parse_options(void);


  //! Setup the available variables
  virtual void do_setup_solution_variables(void);


  //! Solve the MyPoisson equation
  virtual void do_solve(void);


  // ! Print some useful information
  //virtual void do_print_info(void){};

  //! Print some useful information
  virtual void plot_globaldata(void);

  //! We need to create a physical model
  virtual PhysicalModel* create_bulk_model(const ModelOptions& options,
      const Material* mat) const;

  // ! We need to create boundary condition model
//  virtual PhysicalModel* create_boundary_model(const ModelOptions& options,
//      const Material* material_A, const Material* material_B) const{};


  //! We have to provide somehow our solution variables
  virtual void get_solution_secure(const libMesh::Elem* elem,
      std::map<ID, std::vector<double> >& values,
      const std::vector<libMesh::Point>& p);




private:


  //! The constructor
  /*!
   * Being private disables further inheritance.
   */
  Vff(const ModelOptions& options);


  class Options
  {
  public:
    Options(void);

    std::string boundary_conditions;
    std::string substrate_plane;
    double substrate_tol;
    bool substrate_updown;

    double sigma;
    double cutoff;

    //Solver options
    //--------------------------
    //! Minimization method
    std::string method;
    //! Tolerance (eV/A)
    double absolute_tolerance;
    //! Print Level 0, 1, 2 
    int print_lev;
  };

  Options _options;

  Options& get_my_options(void);

  void
  set_boundary(void);

  std::vector<unsigned int> _free_atoms;
  std::vector<unsigned int>& get_free_atoms(void);
  unsigned int _n_free_atoms;
  unsigned int _n_atoms;

// Degree of freedom
  std::vector<double> _dof;
  unsigned int _n_dof;

  // All the coordinates of the atomistic structure, stored in a unique array
  std::vector<double> _coords;
  std::vector<double>& get_coords(void);

  //Initial coordinates. To keep track of hydrogen displacement
  std::vector<double> _initial_coords;

  void
  check_structure(void);

  std::vector<std::vector<double>> _alpha;
  std::vector<std::vector<std::vector<double> > > _beta;
  std::vector<std::vector<double> > _d;
  std::vector<std::vector<std::vector<double> > > _teta;

  void
  resize_parameters(void);

  //! This setup all the arrays of parameters (alpha, beta, d0, teta0)
  void build_parameters(void);

  //! Routine to calculate Keating Potential
  double keating_potential(void);

  //!Routine to calculate analytical gradient of Keating potential
  std::vector<double> keating_gradient(void);


  void optimize(void);

  //! Apply the displacement from temporary coordinates 
  //!to the atomistic structure
  void displace_atoms(void);

  //!Set coordinates used
  void set_coords(void);

  //!Utility for strain projection
  StrainLattice _strain;

  //! Prepare stuff which will be needed for the solutions
  void prepare_solutions(void);

  //! Contains a map of solution for any active element, to speed up the get_element
  std::map<Elem*, Tensor2Gen> _elem_strain_map;
  
 //!Check if we have already calculated the strain tensor, to avoid useless calculation
 bool _has_strain_tensor;

};


inline
Vff::Options&
Vff::get_my_options(void)
{
  return _options;
}

inline
std::vector<unsigned int>&
Vff::get_free_atoms(void)
{
  return _free_atoms;
}

inline
std::vector<double>&
Vff::get_coords(void)
{
  return _coords;
}

inline
const std::vector<double>&
Vff::get_dof(void) const
{
  return _dof;
}

inline
int
Vff::get_n_dof(void) const
{
  return _n_dof;
}
#endif
