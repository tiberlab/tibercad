// $Id$


#ifndef _STRAININTERFACE_H_
#define _STRAININTERFACE_H_

#include "TypeDefs.h"
#include "tiber_dll.h"


#include <set>
#include <vector>
#include <string>

class Elem;
class Point;
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

    
    //! Get the strain in calculation coordinate system
    void get_strain(const Elem* elem, const Point& point, Tensor2Sym& strain);

    //! Get the strain in crystal coordinate system
    void get_crystal_strain(const Elem* elem, const Point& point, Tensor2Sym& strain);

    //! Get the stress in calculation coordinate system
    void get_stress(const Elem* elem, const Point& point, Tensor2Sym& stress);

    //! Get the stress in crystal coordinate system
    void get_crystal_stress(const Elem* elem, const Point& point, Tensor2Sym& stress);


    //! Tells if this interface has a simulation associated
    bool has_simulation(void) const;


    //! Get the associated simulation
    SimulationInterface* get_simulation(void);


  private:

    //! The strain simulation
    SimulationInterface* _simulation;


    //! The ID for variable Strain
    ID _strain_id;

    //! The ID for variable StrainCrystal
    ID _strain_cryst_id;


    //! The ID for variable Stress
    ID _stress_id;

    //! The ID for variable StressCrystal
    ID _stress_cryst_id;

    

    //! Get a tensor
    void _get_data(const Elem* elem, const Point& point, Tensor2Sym& data, ID id);
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
