// $Id$

#ifndef _CONTROL_H_
#define _CONTROL_H_

#include "ModelErrorException.h"
#include "InitFailedException.h"
#include "SolveFailedException.h"

#include <map>

class Device;
class Database;
class SimulationInterface;
class SimulationEnvironment;

//! The control module of TiberCAD
/*!
 * This class is responsible for the setup of the device and
 * the flow control of the simulation
 */
class Control
{

  public:

    //! The constructor
    /*!
     * \param inputfile the input file to be used
     */
    Control(const std::string& inputfile);

    
    //! The destructor
    ~Control(void);

    
    //! Initialize the device and the simulations
    /*!
     * This method calls create_device(), create_materials() and
     * setup_models()
     */
    void init(void) throw (InitFailedException);

    
    //! Runs the simulation
    /*!
     * Calls solve() of the simulations according to the rules given
     * in the input file
     */
    void run_simulation(void) throw (SolveFailedException);
    

    //! Get a reference to the device
    Device& get_device(void);



  private:

    //! Type for the list of simulations
    typedef std::map<const std::string, SimulationInterface*> SimulationMap;

    
    //! Type for the list of simulation environments
    typedef std::map<SimulationInterface*,
            SimulationEnvironment*> EnvironmentMap;

    
    //! The inputfile
    std::string _inputfile;

    
    //! The device we control
    Device* _device;

    
    //! The database we use
    Database* _database;

    
    //! A list of all simulations we control
    SimulationMap _simulations;

    
    //! A list of all simulation environments we control
    EnvironmentMap _simulation_environments;

    
    //! Create the device
    void create_device(void);

    
    //! Create all the materials
    void create_materials(void);


    //! Create and setup the models
    void setup_models(void) throw (ModelErrorException);

};


//
// inline members
//


inline
Control::Control(const std::string& inputfile)
  : _inputfile(inputfile),
    _device(0),
    _database(0)
{
}


inline
Device&
Control::get_device(void)
{
  return *_device;
}



#endif // _CONTROL_H_
