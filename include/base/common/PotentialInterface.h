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

  //! Default constructor creates an empty interface that must be initialized by set_simulation
  PotentialInterface(void);

  //! Constructor that invokes initializations
  PotentialInterface(const std::string& name, const std::string& variable = "ElPotential");


  //! Specify the simulation to use 
  /*!
   * Returns true if \c name refers to a valid potential simulation
   * (examples of potential variables: ElPotential, eQFermi, hQFermi).
   */
  bool set_simulation(const std::string& name, const std::string& variable);


  //! Get the electrostatic potential in specified points
  void get_potential(const libMesh::Elem* elem, const std::vector<libMesh::Point>& p,
      std::vector<double>& potentials, bool local_coord = false);


  //! Get the electrostatic potential in one point
  double get_potential(const libMesh::Elem* elem, const libMesh::Point& p,
      bool local_coord = false);


  //! Get the electrostatic potential at the position of an atom
  double get_potential(const Atom* atom);


  //! Tells if this interface has a simulation associated
  bool has_simulation(void) const;


  //! Get the associated simulation
  SimulationInterface* get_simulation(void);



private:

  //! The potential simulation
  SimulationInterface* _simulation;


  //! The ID as returned from the simulation
  ID _id;

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
