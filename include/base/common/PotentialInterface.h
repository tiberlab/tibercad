// $Id$


#ifndef _POTENTIALINTERFACE_H_
#define _POTENTIALINTERFACE_H_

#include "TypeDefs.h"


#include <set>
#include <vector>
#include <string>

class SimulationInterface;
class Atom;


namespace libMesh
{
  class Elem;
  class Point;
}


//! An interface to a Poisson simulation.
/*!
 * When initialized with a valid simulation name, it will be able
 * to provide a simulated electrostatic potential.
 * The current implementation will return 0 whenever no data can be
 * obtained from the Poisson simulation or if there is no Poisson
 * simulation at all.
 */
class PotentialInterface
{

public:

  //! Default constructor
  PotentialInterface(void);


  //! Specify the Poisson simulation to use
  /*!
   * Returns true if \c name refers to a valid Poisson
   * simulation.
   */
  bool set_simulation(const std::string& name);


  //! Get the electrostatic potential in specified points
  void get_potential(const libMesh::Elem* elem, const std::vector<libMesh::Point>& p,
      std::vector<double>& potentials, bool local_coord = false);


  //! Get the electrostatic potential in one point
  double get_potential(const libMesh::Elem* elem, const libMesh::Point& p,
      bool local_coord = false);


  //! Get the electrostatic potential at the position of an atom
  double get_potential(const Atom* atom);


  //! Get the electron chemical potential in specified points
  void get_el_chem_potential(const libMesh::Elem* elem, const std::vector<libMesh::Point>& p,
      std::vector<double>& potentials, bool local_coord = false);


  //! Get the electron chemical potential in one point
  double get_el_chem_potential(const libMesh::Elem* elem, const libMesh::Point& p,
      bool local_coord = false);


  //! Get the electron chemical potential at the position of an atom
  double get_el_chem_potential(const Atom* atom);


  //! Get the hole chemical potential in specified points
  void get_hl_chem_potential(const libMesh::Elem* elem, const std::vector<libMesh::Point>& p,
      std::vector<double>& potentials, bool local_coord = false);


  //! Get the hole chemical potential in one point
  double get_hl_chem_potential(const libMesh::Elem* elem, const libMesh::Point& p,
      bool local_coord = false);


  //! Get the hole chemical potential at the position of an atom
  double get_hl_chem_potential(const Atom* atom);


  //! Tells if this interface has a simulation associated
  bool has_simulation(void) const;


  //! Get the associated simulation
  SimulationInterface* get_simulation(void);



private:

  //! The potential simulation
  SimulationInterface* _simulation;


  //! The ID as returned from the simulation
  ID _id;

  //! The ID as returned from the simulation
  ID _id_chem_el;

  //! The ID as returned from the simulation
  ID _id_chem_hl;

  //! For efficiency, we will get all potentials together on an atomic position
  //const Atom* _current_atom;

};


//
// inline members
//

inline
bool
PotentialInterface::has_simulation(void) const
{
  return (_simulation == NULL) ? false : true;
}


inline
SimulationInterface*
PotentialInterface::get_simulation(void)
{
  return _simulation;
}


#endif // _POTENTIALINTERFACE_H_
