// $Id$


#ifndef _TEMPERATUREINTERFACE_H_
#define _TEMPERATUREINTERFACE_H_

#include "TypeDefs.h"


#include <set>
#include <vector>
#include <string>

class Elem;
class Point;
class SimulationInterface;


//! An interface to a temperature simulation.
/*!
 * When initialized with a valid simulation name, it will be able
 * to provide a simulated temperature. If there is no simulation,
 * it will just return the temperature from the simulation options
 */
class TemperatureInterface
{

  public:

    //! Default constructor
    TemperatureInterface(void);


    //! Specify the temperature simulation to use
    /*!
     * Returns true if \c name refers to a valid temperature
     * simulation.
     */
    bool set_simulation(const std::string& name);

    
    //! Get the nodal temperatures
    void get_temperature(const Elem* elem, std::vector<double>& temperatures);


    //! Get the temperatures in specified points
    void get_temperature(const Elem* elem, const std::vector<Point>& p,
        std::vector<double>& temperatures);


    //! Get the temperature in one point
    double get_temperature(const Elem* elem, const Point& p);


    //! Tells if this interface has a simulation associated
    bool has_simulation(void) const;


    //! Get the associated simulation
    SimulationInterface* get_simulation(void);



  private:

    //! The temperature simulation
    SimulationInterface* _simulation;


    //! The name of the temperature variable
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
TemperatureInterface::has_simulation(void) const
{
  return (_simulation == NULL) ? false : true;
}


inline
SimulationInterface*
TemperatureInterface::get_simulation(void)
{
  return _simulation;
}


#endif // _TEMPERATUREINTERFACE_H_
