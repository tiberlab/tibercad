// $Id$

#ifndef _SIMULATIONINTERFACE_H_
#define _SIMULATIONINTERFACE_H_

#include "TypeDefs.h"
#include "ModelOptions.h"
#include "InitFailedException.h"
#include "SolveFailedException.h"

// For debugging
#include "reference_counted_object.h"

#include <cassert>
#include <map>
#include <string>

class SimulationEnvironment;
class EquationSystems;

//! The base class for any simulation
class SimulationInterface : public ReferenceCountedObject<SimulationInterface>
{

  public:
    
    //! Destructor
    virtual ~SimulationInterface(void) {};

    //! Set the simulation environment for this simulation
    void set_environment(SimulationEnvironment* environment);

    //! Get the simulation environment for this simulation
    SimulationEnvironment& get_environment(void) const;

    //! Get the ID of this simulation
    ID get_id(void) const;

    //! Get the user defined name of this simulation
    const std::string& get_name(void) const;

    //! Initialize the system
    /*!
     * This method calls do_init() after some health checks
     */
    void init(void) throw (InitFailedException);

    //! Set options for this model
    void set_options(const ModelOptions& options);

    //! Solve the system
    /*!
     * This method calls do_solve() after some health checks
     */
    void solve(void) throw (SolveFailedException);


    //! Create a simulation of type \c type with options
    /*!
     * \param type the simulation type to create
     * \param options the options for the simulation
     * 
     * \return a pointer to the newly created simulation object or \c NULL
     * if the \c typr does not exist.
     */
    static SimulationInterface* create(const std::string& type,
        const ModelOptions& options = ModelOptions());

    //! Destroy a simulation
    /*!
     * \param p the pointer to the simulation to destroy
     */
    static void destroy(SimulationInterface* p);


    //! Find a simulation with name \c name
    /*!
     * \param name the name to look for
     * \return a pointer to the simulation if found, \c NULL otherwise
     */
    static SimulationInterface* find_simulation(const std::string& name);

    //! Check if this simulation is initialized
    bool is_initialized(void) const;

    //! Check if this simulation has already been solved
    /*!
     * This could be useful for models which use results of another
     * simulation.
     */
    bool is_solved(void) const;

    //! Set the relaxation factor
    void set_relaxation_factor(double relax);

    //! Get the relaxation factor
    double get_relaxation_factor(void);


  protected:

    //! Empty constructor
    SimulationInterface(void);

    //! Get a reference to the equation system object
    EquationSystems& get_equation_systems(void) const;

    //! Get the options for this simulator
    ModelOptions& get_options(void);

    //! Do the initialization
    /*!
     * Has to be implemented by derived classes.
     *
     * This method should initialize everything that is needed to do a
     * simulation.
     */
    virtual void do_init(void) = 0;

    //! Do the solve
    /*!
     * Has to be implemented by derived classes.
     *
     * This method does the actual simulation.
     */
    virtual void do_solve(void) = 0;

    //! Parse the options
    /*!
     * This method has to be called \em explicitly somewhere in the derived
     * class. It is not called from \c SimulationInterface::init(), because 
     * in some situations options could change between different calls
     * to \c solve(). It is not called in \c SimulationInterface::solve(),
     * because in other situations this is not necessary. So: call it in
     * \c do_init() or \c do_solve().
     */
    virtual void parse_options(void) = 0;

    //! Get the unique name for the equation system
    const std::string& get_equation_system_name(void) const;


  private:

    //! A typedef for convenience
    typedef std::map<ID, SimulationInterface*> SimulationMap;

    //! The environment for this simulation
    SimulationEnvironment* _environment;

    //! A flag indicating if the simulator is initialized
    bool _is_initialized;

    //! A flag indicating that a simulation has be done
    bool _is_solved;

    //! For self-consistent calculations this could be useful
    double _relaxation_factor;

    //! The ID of this simulation
    /*!
     * The ID is unique for every simulator and is assigned automatically
     * at instantiation.
     */
    ID _id;

    //! Options associated with this model
    /*!
     * These are the options as read from the input file.
     */
    ModelOptions _options;

    //! A user definable name to identify this simulation
    std::string _name;

    //! The unique name for the equation system
    std::string _eq_system_name;

    //! The map containing all simulations with their ID
    static SimulationMap _simulation_map;

    //! create a unique name for the equation system
    void create_equation_system_name(void);

};

inline
void
SimulationInterface::set_environment(SimulationEnvironment* environment)
{
  if (environment != 0)
    _environment = environment;
}

inline
SimulationEnvironment&
SimulationInterface::get_environment(void) const
{
  assert(_environment != 0);

  return *_environment;
}


inline
ID
SimulationInterface::get_id(void) const
{
  return _id;
}


inline
const std::string&
SimulationInterface::get_name(void) const
{
  return _name;
}


inline
const std::string&
SimulationInterface::get_equation_system_name(void) const
{
  return _eq_system_name;
}


inline
void
SimulationInterface::set_options(const ModelOptions& options)
{
  _options += options;
}

inline
ModelOptions&
SimulationInterface::get_options(void)
{
  return _options;
}


inline
bool
SimulationInterface::is_initialized(void) const
{
  return _is_initialized;
}


inline
bool
SimulationInterface::is_solved(void) const
{
  return _is_solved;
}


inline
void
SimulationInterface::set_relaxation_factor(double relax)
{
  _relaxation_factor = relax;
}


inline
double
SimulationInterface::get_relaxation_factor(void)
{
  return _relaxation_factor;
}


#endif // _SIMULATIONINTERFACE_H_
