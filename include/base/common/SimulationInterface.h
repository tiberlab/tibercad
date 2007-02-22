// $Id$

#ifndef _SIMULATIONINTERFACE_H_
#define _SIMULATIONINTERFACE_H_

#include "TypeDefs.h"
#include "ModelOptions.h"
#include "InitFailedException.h"
#include "SolveFailedException.h"
#include "ModelErrorException.h"

// LibMesh includes
// For debugging
#include "reference_counted_object.h"
#include "numeric_vector.h"

#include <cassert>
#include <map>
#include <set>
#include <string>

class SimulationEnvironment;
class EquationSystems;
class PhysicalModel;
class BoundaryProperties;
class Control;
class Mesh;


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


    //! Get the default name for this simulation
    /*!
     * This can be useful to identify a named simulation for which
     * we dont't know the name neither the type
     */
    std::string get_default_name(void) const;


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


    //! Create a physical model that can be used by this type of simulation
    /*!
     * The default behaviour defined in the base class is to return the
     * NULL pointer, because there could be simulations that don't need
     * a physical model. If a derived class reimplements this method (which
     * will normally be the case) it should notify about errors by throwing
     * a ModelErrorException.
     */
    virtual PhysicalModel*
      create_physical_model(const ModelOptions& options) const
      throw (ModelErrorException);

    
    //! Create a boundary model that can be used by this type of simulation
    /*!
     * The default behaviour defined in the base class is to return the
     * NULL pointer, because there could be simulations that don't need
     * a boundary model. If a derived class reimplements this method (which
     * will normally be the case) it should notify about errors by throwing
     * a ModelErrorException.
     */
    virtual BoundaryProperties*
      create_boundary_model(const ModelOptions& options) const
      throw (ModelErrorException);


    //! Initialize the system
    /*!
     * This method calls do_init() after some health checks
     */
    void init(void) throw (InitFailedException);

    
    //! Set options for this model
    /*!
     * \param options the options to be changed
     *
     * This method will take all options given in \c options and update
     * the options of the solver accordingly. It does not touch options
     * which are not present in the argument of this method
     */
    void set_options(const ModelOptions& options);

    
    //! Solve the system
    /*!
     * This method calls do_solve() after some health checks
     */
    void solve(void) throw (SolveFailedException);


    //! Write results to file
    void plot(void);


    //! Find a simulation with name \c name
    /*!
     * \param name the name to look for
     * \return a pointer to the simulation if found, \c NULL otherwise
     *
     * \c name can be one of the following:
     * \li the user defined name of a simulation
     * \li the identifier of the simulation as used for creation
     * \li the empty string
     *
     * In the second case, the first simulation of this type will be
     * returned. In the third case, the first simulation will be returned.
     * 
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
    double get_relaxation_factor(void) const;

    
    //! Get a pointer to the solution vector
    NumericVector<Real>& get_solution_vector(void);

    
    /*!
     * \copydoc build_nodal_results()
     *
     * Calls build_nodal_results()
     */
    void get_nodal_results(const std::set<std::string>& variables,
        std::vector<double>& results, std::vector<std::string>& legend);


    /*!
     * \copydoc build_elemental_results()
     *
     * Calls build_elemental_results()
     */
    void get_elemental_results(const std::set<std::string>& variables,
        std::vector<double>& results, std::vector<std::string>& legend);


    //! Build a vector with some integrated quantities
    /*!
     * Calls build_integrated_quantities()
     *
     * This method is used in sweeps for creating files like I-V
     * characteristics and similar.
     *
     * \param names the identifier of the quantities to plot
     * \param values the vector where the values will be put
     */
    void get_integrated_quantities(const std::set<std::string>& names,
        std::vector<double>& values);
    

    //! Get the description for some integrated quantities
    /*!
     * calls build_integrated_quantities_description()
     *
     * cf. get_integrated_quantities()
     *
     * \param names the identifier of the quantities to plot
     * \param legend the legend for the plot values, has usually the same
     * size as \c values in get_integrated_quantities()
     * \param description a description for each of the known quantities
     */
    void get_integrated_quantities_description(
        const std::set<std::string>& names,
        std::vector<std::string>& legend,
        std::vector<std::string>& description);


    //! Get the type of this simulation
    /*!
     * The type is the identifying string which defines at creation time
     * which simulation to create. It's the same string one writes in the
     * input file.
     */
    const std::string& get_type(void) const;


    //! Set the Control object
    void set_control(Control* control);


    //! Get the Control object
    Control& get_control(void);



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
    

    //! Build nodal result vector for the given variables
    /*!
     * \param variables the identifier for the quantities that should be
     * putted into the vector
     * \param results the vector that will contain the results
     * \param legend the legend for the values in \c results
     */
    virtual void build_nodal_results(const std::set<std::string>& variables,
        std::vector<double>& results, std::vector<std::string>& legend);


    //! Build elemental result vector for the given variables
    /*!
     * \param variables the identifier for the quantities that should be
     * putted into the vector
     * \param results the vector that will contain the results
     * \param legend the legend for the values in \c results
     */
    virtual void build_elemental_results(const std::set<std::string>& variables,
        std::vector<double>& results, std::vector<std::string>& legend);


    //! Build a vector with some integrated quantities
    /*!
     * \param names the identifier for the quantities that should be
     * putted into the vector
     * \param values the vector that will contain the results
     * \param legend the legend for the values in \c results
     *
     * This method is thought for quantities that do not depend on coordinates
     * but are e.g. scalar integrated quantities like terminal currents.
     * They can also be scalar values depending on some other quantity like
     * energy (e.g. for some spectrum). In this case, the legend should contain
     * the corresponding values of the independent variable.
     */
    virtual void build_integrated_quantities(
        const std::set<std::string>& names,
        std::vector<double>& values) {};

    //! Create legend and description for integrated quantities
    /*!
     * cf. build_integrated_quantities()
     *
     * The return values of this method are used in printing data files
     */
    virtual void build_integrated_quantities_description(
        const std::set<std::string>& names,
        std::vector<std::string>& legend,
        std::vector<std::string>& description) {};



  private:

    //! A typedef for convenience
    typedef std::map<ID, SimulationInterface*> SimulationMap;

    
    //! The environment for this simulation
    SimulationEnvironment* _environment;


    //! The Control object which controls this simulation
    Control* _control;

    
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


    //! The identifying string for the type of this simulation
    std::string _type;

    
    //! The unique name for the equation system
    std::string _eq_system_name;

    
    //! The map containing all simulations with their ID
    static SimulationMap _simulation_map;

    
    //! create a unique name for the equation system
    void create_equation_system_name(void);


    //! Set the simulation type (= identifier)
    /*!
     * The identifier is used at creation time to know which type of
     * simulation to create.
     */
    void set_type(const std::string& type);


};


//
// inline members
// 


inline
void
SimulationInterface::set_control(Control* control)
{
  _control = control;
}


inline
Control&
SimulationInterface::get_control(void)
{
  assert(_control != NULL);
  return *_control;
}



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
SimulationInterface::get_relaxation_factor(void) const
{
  return _relaxation_factor;
}



inline
void
SimulationInterface::set_type(const std::string& type)
{
  _type = type;
}



inline
const std::string&
SimulationInterface::get_type(void) const
{
  return _type;
}


inline
void
SimulationInterface::get_elemental_results(
    const std::set<std::string>& variables,
    std::vector<double>& results, std::vector<std::string>& legend)
{
  results.resize(0);
  legend.resize(0);
  build_elemental_results(variables, results, legend);
}


inline
void
SimulationInterface::get_nodal_results(
    const std::set<std::string>& variables,
    std::vector<double>& results, std::vector<std::string>& legend)
{
  results.resize(0);
  legend.resize(0);
  build_nodal_results(variables, results, legend);
}



inline
void
SimulationInterface::get_integrated_quantities(
    const std::set<std::string>& names, std::vector<double>& values)
{
  values.resize(0);
  build_integrated_quantities(names, values);
}


inline
void
SimulationInterface::get_integrated_quantities_description(
    const std::set<std::string>& names,
    std::vector<std::string>& legend,
    std::vector<std::string>& description)
{
  legend.resize(0);
  description.resize(0);
  build_integrated_quantities_description(names, legend, description);
}



inline
void
SimulationInterface::build_nodal_results(
    const std::set<std::string>& variables,
    std::vector<double>& results, std::vector<std::string>& legend)
{
  ignore_unused_variable(variables);
  ignore_unused_variable(results);
  ignore_unused_variable(legend);
}



inline
void
SimulationInterface::build_elemental_results(
    const std::set<std::string>& variables,
    std::vector<double>& results, std::vector<std::string>& legend)
{
  ignore_unused_variable(variables);
  ignore_unused_variable(results);
  ignore_unused_variable(legend);
}


#endif // _SIMULATIONINTERFACE_H_
