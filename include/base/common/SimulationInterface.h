// $Id$

#ifndef _SIMULATIONINTERFACE_H_
#define _SIMULATIONINTERFACE_H_

#include "tiber_config.h"
#include "TypeDefs.h"
#include "ModelOptions.h"
#include "InitFailedException.h"
#include "SolveFailedException.h"
#include "ModelErrorException.h"
#include "Scaling.h"
#include "FiniteElement.h"

// LibMesh includes
// For debugging
#include "reference_counted_object.h"
#include "numeric_vector.h"

#include <cassert>
#include <map>
#include <set>
#include <string>


#ifdef BUILD_TIBER_MODULES
/*!
 * \def TIBER_MODULE(classname, libname)
 *
 * \brief Creates methods to create and destroy a simulation object
 * 
 * In each implementation derived from SimulationInterface, put
 * this macro somewhere in the source file to be able to compile
 * it as TiberCad module.
 *
 * \param name the name of the class that should be 'creatable'
 * \param libname the name for this module
 *
 * \c libname will be used to create the library name, and the model
 * will have to be referred to in the input file by \c libname
 */
# ifndef TIBER_MODULE
#  define TIBER_MODULE(classname, libname) \
  extern "C" { \
    void destroy(SimulationInterface* p) { \
      delete p; \
    } \
    classname* create(void) { \
      return new classname(); \
    } \
    const char* _tiber_module_ ## libname = #libname; \
    const char* library_name(void) { \
      return _tiber_module_ ## libname; \
    } \
  }
# endif
#else
# ifndef TIBER_MODULE
#  define TIBER_MODULE(classname, libname)
# endif
#endif


class SimulationEnvironment;
class EquationSystems;
class PhysicalModel;
class BoundaryProperties;
class Control;
class Mesh;
class Point;

//! The base class for any simulation
class SimulationInterface : public ReferenceCountedObject<SimulationInterface>
{

  public:
    
    //! Destructor
    virtual ~SimulationInterface(void);


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


    //! Solve the system for equilibrium
    /*!
     * Calls do_equilibrium()
     *
     */
    void solve_equilibrium(void) throw (SolveFailedException);
    

    //! Solve the system
    /*!
     * This method calls do_solve() after some health checks
     */
    void solve(void) throw (SolveFailedException);


    //! Write results to file
    /*!
     * This method calls do_plot() after some health checks
     */
    void plot(void);


    //! Remember current solution (on the current mesh!)
    /*!
     * Calls do_remember_current_solution()
     */
    /*!
     * If \c id is invalid or not specified, a new ID will be generated and
     * the solution stored with this new ID. Otherwise a solution stored at 
     * \c id will be overwritten.
     */
    ID remember_current_solution(ID id = 0);


    //! Set to the remembered solution number \c id
    /*!
     * Calls do_set_to_remembered_solution()
     */
    void set_to_remembered_solution(ID id);


    //! Clear a remembered solution
    /*!
     * Calls do_delete_remembered_solution()
     */
    void delete_remembered_solution(ID id);


    //! Get a remembered solution
    NumericVector<double>* get_remembered_solution(ID id);


    /*!
     * \brief Get the maximum norm of the difference between the
     * current and a remembered solution.
     *
     * Calls do_maximum_norm_of_difference()
     *
     */
    double get_maximum_norm_of_difference(ID id);


    //! Scale the current solution by a real factor
    /*!
     * Calls do_scale_solution()
     */
    void scale_solution(double factor);


    //! Adds a scaled remembered solution to the current solution
    /*!
     * Calls do_add_scaled_remembered_solution()
     */
    void add_scaled_remembered_solution(ID id, double factor);


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


    //! Get a pointer to the simulation with ID \c id
    static SimulationInterface* get_simulation(ID id);

    
    //! Check if this simulation is initialized
    bool is_initialized(void) const;

    
    //! Check if this simulation has already been solved
    /*!
     * This could be useful for models which use results of another
     * simulation.
     */
    bool is_solved(void) const;


    //! Check if the equilibrium for this simulation has been calculated
    bool equilibrium_done(void) const;
    
    
    //! Set the relaxation factor
    //void set_relaxation_factor(double relax);

    
    //! Get the relaxation factor
    //double get_relaxation_factor(void) const;




    /*!
     * \copydoc convert_variable_name_to_id()
     *
     */
    ID get_variable_id(const std::string& variable_name) const;


    //! Get solution values on the nodes of a specified element
    /*!
     *
     * \param elem a pointer to the element
     * \param id identifier for the variable to be returned
     * \param values a vector to store the values. The vector index corresponds to
     * the node.
     *
     * \return false if no data can be found for \c elem
     */
    bool get_solution(const Elem* elem, ID id,
        std::vector<double>& values);


    //! Get solution values on the nodes of a specified element
    /*!
     *
     * \param elem a pointer to the element
     * \param ids identifiers for the variables to be returned
     * \param values a vector to store the values. The vector index
     * corresponds to the node.
     *
     * \return false if no data can be found for \c elem
     */
    bool get_solution(const Elem* elem, const std::set<ID>& ids,
        std::vector<std::map<ID, double> >& values);


    //! Get solution values on one inner point of a specified element
    /*!
     *
     * \param elem a pointer to the element
     * \param p the point (assumed to lie in \c elem)
     * \param ids identifiers for the variables to be returned
     * \param values a vector to store the values. The vector index
     * corresponds to the point.
     *
     * \return false if no data can be found for \c elem
     */
    bool get_solution(const Elem* elem, const Point& p,
        const std::set<ID>& ids, std::map<ID, double>& values);
  

    //! Get solution values on one inner point of a specified element
    /*!
     *
     * \param elem a pointer to the element
     * \param p the point (assumed to lie in \c elem)
     * \param id identifier for the variable to be returned
     * \param value the output value
     *
     * \return false if no data can be found for \c elem
     */
    bool get_solution(const Elem* elem, const Point& p, ID id, double& value);
  

    //! Get solution values on inner points of a specified element
    /*!
     *
     * \param elem a pointer to the element
     * \param p a vector with the points. All Points are assumed to lie in
     * \c elem
     * \param id identifier for the variable to be returned
     * \param values a vector to store the values. The vector index
     * corresponds to the point.
     *
     * \return false if no data can be found for \c elem
     */
    bool get_solution(const Elem* elem, const std::vector<Point>& p,
        ID id, std::vector<double>& values);
  

    //! Get solution values on inner points of a specified element
    /*!
     *
     * \param elem a pointer to the element
     * \param p a vector with the points. All Points are assumed to lie in
     * \c elem
     * \param ids identifiers for the variables to be returned
     * \param values a vector to store the values. The vector index
     * corresponds to the point.
     *
     * \return false if no data can be found for \c elem
     */
    bool get_solution(const Elem* elem, const std::vector<Point>& p,
        const std::set<ID>& ids, std::vector<std::map<ID, double> >& values);
    



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
    void get_integrated_quantities(const std::set<std::string>& variables,
        std::vector<double>& values);
    

    //! Get the description for some integrated quantities
    /*!
     * calls build_integrated_quantities_description()
     *
     * cf. get_integrated_quantities()
     *
     * \param variables the identifier of the quantities to plot
     * \param legend the legend for the plot values, has usually the same
     * size as \c values in get_integrated_quantities()
     * \param description a description for each of the known quantities
     *
     * \note This method must not access the mesh or DOF maps or similar
     * as it is not guaranteed the the mesh is prepared for the currently
     * accessed simulation.
     */
    void get_integrated_quantities_description(
        const std::set<std::string>& variables,
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


    //! Set the scaling parameters
    void set_scaling(const Scaling& scaling);


    //! Get the scaling parameters
    Scaling& get_scaling(void);


    //! Build a finite element
    /*!
     * Currently only Lagrange elements are supported.
     *
     * \param dim the space dimension for this element
     * \param type the type of finite element
     * \param scale include length scaling if true
     */
    AutoPtr<FEBase> build_finite_element(unsigned int dim, FEType type,
        bool scale = false);



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


    //! Solve for equilibrium
    /*!
     * Can be implemented in derived classes to solve for some
     * equilibrium state. This method will be called for any simulation
     * before calling solve().
     *
     * Use it e.g. to solve the Poisson equation for equilibrium, to set
     * some reasonable starting values etc.
     */
    virtual void do_equilibrium(void) {};

    
    //! Do the solve
    /*!
     * Has to be implemented by derived classes.
     *
     * This method does the actual simulation.
     */
    virtual void do_solve(void) = 0;


    //! Do a plot
    /*!
     * This method creates files with the results of the simulation.
     * In almost all cases the default implementation should be ok, but
     * in some it is not, especially in "compound" simulations as
     * sweeps, selfconsitent solvers etc.
     */
    virtual void do_plot(void);

    
    //! Remember the current solution
    /*!
     * The default action is to clone the solution vector of the 
     * equations system. If \c id is invalid, a new ID will be generated
     * and returned, if not, the remembered solution \c id will be overwritten.
     * If this method is reimplemented in a derived class, it should behave in
     * the same way!
     */
    virtual ID do_remember_current_solution(ID id = 0);


    //! Return true if the system has a solution vector
    bool has_solution_vector(void);

    
    //! Get a pointer to the solution vector
    virtual NumericVector<double>& get_solution_vector(void);


    //! Set to the remembered solution number \c id
    virtual void do_set_to_remembered_solution(ID id);


    //! Delete a remembered solution
    virtual void do_delete_remembered_solution(ID id);
    

    //! Build the maximum norm of the solution difference
    /*!
     * This method returns the maximum norm \f$\Vert x - x_i\Vert_\infty\f$
     * of the difference between the current solution \f$x\f$ and some
     * remembered solution \f$x_i\f$
     *
     * \param id the id of the remembered solution.
     *
     * If the remembered solution doesn't exist, return value is zero.
     *
     * The default action should be ok in most cases. But it could be useful
     * to reimplement this method when a simulation uses some scaling
     */
    virtual double do_maximum_norm_of_difference(ID id);


    //! Scale the current solution
    /*!
     * \param factor the scaling factor
     */
    virtual void do_scale_solution(double factor);
    

    //! Add a scaled remembered solution
    /*!
     * \param id the id of the remembered solution
     * \param factor the scaling factor
     */
    virtual void do_add_scaled_remembered_solution(ID id, double factor);


    //! Parse the options
    /*!
     * This method has to be called \em explicitly somewhere in the derived
     * class. It is not called from \c SimulationInterface::init(), because 
     * in some situations options could change between different calls
     * to \c solve(). It is \em not called in \c SimulationInterface::solve(),
     * because in other situations this is not necessary. So: call it in
     * \c do_init() or \c do_solve().
     */
    virtual void parse_options(void) = 0;

    
    //! Get the unique name for the equation system
    const std::string& get_equation_system_name(void) const;

    
    //! Get the ID for a given variable name
    /*!
     * An ID of \c INVALID_ID represents an unknown variable.
     *
     * \param the name of the variable as a string
     * \return the ID of the variable as a numerical value
     */
    virtual ID convert_variable_name_to_id(const std::string& variable_name) const;
    
    
    //! Get solution values on the nodes of a specified element
    /*!
     *
     * \param elem a pointer to the element
     * \param ids identifiers for the variables to be returned
     * \param values a vector to store the values. The first index
     * corresponds to the node.
     *
     * \pre 
     * \li \c elem is an active element of this simulation
     * \li values is already resized to the number of nodes
     */
    virtual void get_solution_secure(const Elem* elem,
        const std::set<ID>& ids, std::vector<std::map<ID, double> >& values) {};


    //! Get solution values on inner points of a specified element
    /*!
     *
     * \param elem a pointer to the element
     * \param p a vector with the points. All Points are assumed to lie in
     * \c elem
     * \param ids identifiers for the variables to be returned
     * \param values a vector to store the values. The first index
     * corresponds to the point.
     * 
     * \pre 
     * \li \c elem is an active element of this simulation
     * \li \c elem contains all points of \c p
     * \li values is already resized to the number of points
     */
    virtual void get_solution_secure(const Elem* elem,
        const std::vector<Point>& p, const std::set<ID>& ids,
        std::vector<std::map<ID, double> >& values) {};

    
    //! Get the solution values on one point of an element
    /*!
     *
     * \param elem a pointer to the element
     * \param p the point (assumed to lie in \c elem)
     * \param values a vector to store the values.
     * 
     * \pre 
     * \li \c elem is an active element of this simulation
     * \li \c elem contains all points of \c p
     * \li values is already resized to the number of points
     */
    virtual void get_solution_secure(const Elem* elem,
        const Point& p, const std::set<ID>& ids,
        std::map<ID, double>& values);




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
     * \param variables the identifier for the quantities that should be
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
        const std::set<std::string>& variables,
        std::vector<double>& values) {};

    //! Create legend and description for integrated quantities
    /*!
     * cf. build_integrated_quantities()
     *
     * The return values of this method are used in printing data files
     */
    virtual void build_integrated_quantities_description(
        const std::set<std::string>& variables,
        std::vector<std::string>& legend,
        std::vector<std::string>& description) {};



  private:

    //! A typedef for convenience
    typedef std::map<ID, SimulationInterface*> SimulationMap;
    

    //! The creation method signature
    typedef SimulationInterface* (*create_t)(void);

    
    //! The destruction method signature
    typedef void (*destroy_t)(SimulationInterface*);


    //! The type for library handles
    typedef void* libhandle_t;


    
    //! The environment for this simulation
    SimulationEnvironment* _environment;


    //! The Control object which controls this simulation
    Control* _control;

    
    //! A flag indicating if the simulator is initialized
    bool _is_initialized;

    
    //! A flag indicating that a simulation has been done
    bool _is_solved;

    
    //! A flag indicating that equilibrium has been done
    bool _equilibrium_is_solved;


    //! For self-consistent calculations this could be useful
    //double _relaxation_factor;



    //! The library handle for this simulation type
    libhandle_t _libhandle;

    
    //! The creation method for this simulation type
    create_t _create;

    
    //! The destruction method for this simulation type
    destroy_t _destroy;


    
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


    //! The scaling parameters
    Scaling _scaling;

    
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


    //! A map with remembered solutions
    std::map<ID, NumericVector<double>*> _remembered_solutions;



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
SimulationInterface*
SimulationInterface::get_simulation(ID id)
{
  SimulationInterface* sim = NULL;
  
  SimulationMap::iterator it(_simulation_map.find(id));
  if (it != _simulation_map.end())
    sim = it->second;

  return sim;
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


/* inline */
/* const ModelOptions& */
/* SimulationInterface::see_options(void) */
/* { */
/*   return _options; */
/* } */


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
bool
SimulationInterface::equilibrium_done(void) const
{
  return _equilibrium_is_solved;
}


/*
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
*/


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


inline
ID
SimulationInterface::remember_current_solution(ID id)
{
  return do_remember_current_solution(id);
}


inline
void
SimulationInterface::set_to_remembered_solution(ID id)
{
  do_set_to_remembered_solution(id);
}


inline
void
SimulationInterface::delete_remembered_solution(ID id)
{
  do_delete_remembered_solution(id);
}


inline
double
SimulationInterface::get_maximum_norm_of_difference(ID id)
{
  return do_maximum_norm_of_difference(id);
}


inline
void
SimulationInterface::scale_solution(double factor)
{
  do_scale_solution(factor);
}


inline
void
SimulationInterface::add_scaled_remembered_solution(ID id, double factor)
{
  do_add_scaled_remembered_solution(id, factor);
}



inline
Scaling&
SimulationInterface::get_scaling(void)
{
  return _scaling;
}


inline
void
SimulationInterface::set_scaling(const Scaling& scaling)
{
  _scaling = scaling;
}


inline
ID
SimulationInterface::get_variable_id(const std::string& variable_name) const
{
  return convert_variable_name_to_id(variable_name);
}


inline
ID
SimulationInterface::convert_variable_name_to_id(
    const std::string& variable_name) const
{
  ignore_unused_variable(variable_name);
  return INVALID_ID;
}


inline
void
SimulationInterface::get_solution_secure(const Elem* elem,
    const Point& p, const std::set<ID>& ids, std::map<ID, double>& values)
{
  std::vector<Point> pvec(1, p);
  std::vector<std::map<ID, double> > valvec(1);
  
  get_solution_secure(elem, pvec, ids, valvec);

  values = valvec[0];
}




#endif // _SIMULATIONINTERFACE_H_
