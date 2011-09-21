#ifndef _TIGHTBINDING_H_
#define _TIGHTBINDING_H_

//-----------------------------------------------------------------------------------------

#include "EigenvalueProblem.h"
#include "SimulationEnvironment.h"
#include "Database.h"
#include "Specie.h"


//forward declaration
class Device;
class DftbpWrapper;
class UptWrapper;

//!Main class for Atomistic Tight Binding simulation at equilibrium
//! DFTB code is used for simulations
class TightBinding : public EigenvalueProblem{


public:

  enum Unused
  {
    UNKNOWN = 0,
    CHARGE
  };

  enum Shell
  {
    NONE = 0,
    S = 1,
    P = 2,
    D = 3
  };


  //! Constructor
  TightBinding(const ModelOptions& options);

  //! Destructor
  ~TightBinding();

  //! Create TightBinding object
  static TightBinding* create(const ModelOptions& options);

  virtual BoundaryProperties* create_boundary_model(const ModelOptions &options) const
  throw (ModelErrorException);

  virtual void get_solution_secure(const Elem* elem,
      const std::set<ID>& ids, std::vector<std::map<ID, double> >& values);

  virtual void
  get_solution_secure(const Elem* elem, const std::vector<Point>& p,
      const std::set<ID>& ids, std::vector<std::map<ID, double> >& values);
 
//  //! Order the solution in correct mode
//  virtual void build_elemental_results(const std::set<std::string>& variables,
//      std::vector<double>& results, std::vector<std::string>& legend);

  

private:


 
protected:

  virtual void do_init (void);

  virtual void do_solve (void);

  virtual void parse_options(void);

  //! Get hubbard parameters for species composing atomistic structure
  /*!
   *  (it fills a class member map _u_hub)
   */
  virtual void obtain_hubbard_parameters(void);


  //! Map of map containing hubbard parameters for any specie and any shell
  /*!
   * Usage: _u_hub[<specie>][shell] = hubbard_index
   */
  std::map<Specie, std::map<Shell, double> > _u_hub;

  //! Build charge density on given point
  double build_rho(const Point& r);

  // Build states density 
//  virtual void 
  //build_statedens(std::vector<double>& values, const Point& r){};

  //! Charge variation (Mulliken Analisys) on each atom
  std::vector<double> _mulliken_netcharges;

  //! Build a vector of potential projection over atom orbitals
  void project_potential(const std::string providing_model, const std::string mode);

  //!Vector for atom-projected potential shifts
  std::vector<double> _pot_shift;

  //! minimum potential
  double _pot_min;

  //!Vector for atom-projected electron chemical potential
  std::vector<double> _el_chem_pot;

  //!Vector for atom-projected hole chemical potential
  std::vector<double> _hl_chem_pot;

  std::vector< std::vector<unsigned int> > _elem_to_atoms;

  void build_map_elem_atoms(double projection_length);

  const std::vector<unsigned int>& get_neigh_atoms(unsigned int id);

};


inline
TightBinding* TightBinding::create(const ModelOptions& options)
{
  return new  TightBinding(options);
}

inline
const std::vector<unsigned int>& TightBinding::get_neigh_atoms(unsigned int id)
{
  return _elem_to_atoms[id];
}


#endif
