// $Id$


#ifndef _POTENTIALINTERFACE_H_
#define _POTENTIALINTERFACE_H_

#include "TypeDefs.h"


#include <set>
#include <vector>
#include <string>

class Elem;
class Point;
class SimulationInterface;


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


  //! Get the nodal electrostatic potentials
  void get_potential(const Elem* elem, std::vector<double>& potentials);


  //! Get the electrostatic potential in specified points
  void get_potential(const Elem* elem, const std::vector<Point>& p,
      std::vector<double>& potentials);


  //! Get the electrostatic potential in one point
  double get_potential(const Elem* elem, const Point& p);


  //THIS SECTION IS TEMPORARY FOR IWCE CALCULATION
  //------------------------------------------------------------
  //! Get the nodal electron chemical potentials
  void get_el_chem_potential(const Elem* elem, std::vector<double>& potentials);


  //! Get the electron chemical potential in specified points
  void get_el_chem_potential(const Elem* elem, const std::vector<Point>& p,
      std::vector<double>& potentials);


  //! Get the electron chemical potential in one point
  double get_el_chem_potential(const Elem* elem, const Point& p);

  //! Get the nodal hole chemical potentials
  void get_hl_chem_potential(const Elem* elem, std::vector<double>& potentials);


  //! Get the hole chemical potential in specified points
  void get_hl_chem_potential(const Elem* elem, const std::vector<Point>& p,
      std::vector<double>& potentials);


  //! Get the hole chemical potential in one point
  double get_hl_chem_potential(const Elem* elem, const Point& p);

  //! The ID as returned from the simulation
  ID _id_chem_el;

  //! The ID as returned from the simulation
  ID _id_chem_hl;
  //-------------------------------------------------------------


  //! Tells if this interface has a simulation associated
  bool has_simulation(void) const;


  //! Get the associated simulation
  SimulationInterface* get_simulation(void);



private:

  //! The potential simulation
  SimulationInterface* _simulation;


  //! The name of the potential variable
  static std::string _variable_name;


  //! The ID as returned from the simulation
  ID _id;


  //! We need the ID in a set for the function calls
  std::set<ID> _id_set;

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
