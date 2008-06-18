// $Id$


#ifndef _STRAININTERFACE_H_
#define _STRAININTERFACE_H_

#include "TypeDefs.h"

#include "vector_value.h"

#include <set>
#include <vector>
#include <string>

class Elem;
class SimulationInterface;
class Tensor2Sym;


//! An interface to a strain simulation.
/*!
 * When initialized with a valid simulation name, it will be able
 * to provide a simulated strain. If there is no simulation,
 * it will just return 0 strain
 */
class StrainInterface
{

  public:

    //! Default constructor
    StrainInterface(void);


    //! Specify the strain simulation to use
    /*!
     * Returns true if name corresponds to an existing strain
     * simulation.
     */
    bool set_simulation(const std::string& name);

    
    //! Get the strain on an element
    void get_strain_data(const Elem* elem, Tensor2Sym& strain,
        RealVectorValue& polarization);


    //! Tells if this interface has a simulation associated
    bool has_simulation(void) const;


    //! Get the associated simulation
    SimulationInterface* get_simulation(void);


  private:

    //! The strain simulation
    SimulationInterface* _simulation;


    //! The ID as returned from the simulation
    ID _id;
    std::vector<ID> _strain_ids;

    
    //! We need the ID in a set for the function calls
    std::set<ID> _id_set;

};


//
// inline members
//

inline
bool
StrainInterface::has_simulation(void) const
{
  return (_simulation == NULL) ? false : true;
}


inline
SimulationInterface*
StrainInterface::get_simulation(void)
{
  return _simulation;
}


#endif // _STRAININTERFACE_H_
