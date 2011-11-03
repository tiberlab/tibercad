// $Id$

/*! \file SimulationInterface.h */

#ifndef _SIMULATIONINTERFACE_H_
#define _SIMULATIONINTERFACE_H_

#include "tiber_config.h"
#include "TiberModelObject.h"
#include "TypeDefs.h"
#include "HashMap.h"
#include "IDSet.h"
#include "SolutionDescriptor.h"
#include "InitFailedException.h"
#include "SolveFailedException.h"
#include "ModelErrorException.h"
#include "RuntimeException.h"
#include "Scaling.h"
#include "FiniteElement.h"

// LibMesh includes
#include "numeric_vector.h"

#include <cassert>
#include <map>
#include <set>
#include <string>

#include "TiberModule.h"

#ifndef TB_MAX_SIM
#define TB_MAX_SIM 100
#endif

class SimulationEnvironment;
class Embracing;
class Device;
class EquationSystems;
class TiberEqSystem;
class PhysicalModel;
class BoundaryProperties;
class Material;
class MaterialBoundary;
class EdgeObject;
class NodeObject;
class MeshBase;
class Point;
class DofObject;
class Atom;
class AtomisticStructure;


//! The base class for any simulation
class SimulationInterface : public TiberModelObject
{

  private:

    //! A typedef for convenience
    typedef std::map<ID, SimulationInterface*> SimulationMap;


  public:


    //! An iterator to iterate over all simulations
    class SimulationIterator
    {
      public:
        SimulationIterator(void) : _iter(_simulation_map.end()) {}
        SimulationIterator(const SimulationIterator& it) : _iter(it._iter) {}
        SimulationIterator(const SimulationMap::iterator& it) : _iter(it) {}

        SimulationIterator& operator++(void)
        {
          ++_iter;
          return *this;
        }

        SimulationIterator& operator=(const SimulationIterator& it)
        {
          _iter = it._iter;
          return *this;
        }

        bool operator==(const SimulationIterator& it)
        {
          return (_iter == it._iter);
        }

        bool operator!=(const SimulationIterator& it)
        {
          return !this->operator==(it);
        }

        SimulationInterface* operator*(void)
        {
          return _iter->second;
        }


      private:
        SimulationMap::iterator _iter;
    };



    //! Destructor
    virtual ~SimulationInterface(void);


    //! Setup the environment for this simulation
    /*!
     * The SimulationEnvironment for each simulation is unique
     */
    void setup_environment(Device& device, const std::set<ID>& region_numbers);


    //! Check if the simulation has a valid environment
    bool has_environment(void) const;


    //! Get the simulation environment for this simulation
    SimulationEnvironment& get_environment(void) const;


    //! Get the ID of this simulation
    ID get_id(void) const;


    //! Returns true if this is a task module
    /*!
     * A task module is a simulation module which performs
     * tasks like sweep, selfconsistent cycle etc. on other
     * simulation modules. Such a module typically does not
     * have an associated environment, and may have no
     * solution variables.
     *
     * A module that is a task module has to declare this
     * explicitly by calling \c is_task(true).
     */
    bool is_task(void) const;


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


    //! Destroy a SimulationInterface object
    /*!
     * Checks if the given pointer is valid
     */
    static void destroy(SimulationInterface* sim);


    //! Get the iterator to the first simulation
    static SimulationIterator simulations_begin(void);


    //! Get the past-the-end iterator of the simulation list
    static SimulationIterator simulations_end(void);


    //! Prepare the simulation
    /*!
     * This does some preliminary initialization of data that might be needed
     * during initialization of models or other modules.
     */
    void prepare(void);


    //! Initialize the system
    /*!
     * This method calls do_init() after some health checks
     * and print_info();
     */
    void init(void);


    //! Reinitialize the system before each solve
    /*!
     * This method is called before calling do_solve()
     */
    void reinit(void);


    //! Setup the available solution variables
    /*!
     * This has to be done before calling init() of all models
     * because others could want to know about solution variables
     * during initialisation.
     */
    void setup_solution_variables(void);


    //! Print some useful info
    /*!
     * Calls do_print_info()
     */
    void print_info(void);


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
    void solve(void);


    //! Obtain the solve sequence number
    /*!
     * Every call to solve increments the solve sequence number
     * by one. This number can be used to decide if operations
     * that depend on results of this model need to be carried out.
     */
    unsigned int get_solve_sequence_number(void) const;


    //! Write results to file
    /*!
     * This method calls do_plot() after some health checks
     */
    void plot(void);


    //! Save the current state
    void save_state(const std::string& file = "");


    //! Load a saved state
    /*!
     * \return true if a state has been loaded
     */
    bool load_state(const std::string& file = "");


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

    /*!
     * \brief Get the l2 norm of the difference between the
     * current and a remembered solution.
     *
     * Calls do_l2_norm_of_difference()
     *
     */
    double get_l2_norm_of_difference(ID id);



    //! Return true if the system has a solution vector
    bool has_solution_vector(void);


    //! Get a pointer to the solution vector
    /*!
     * Calls do_get_solution_vector()
     */
    NumericVector<double>& get_solution_vector(void);


    //! Set the current solution to a given value
    /*!
     * Calls do_set_solution_vector()
     */
    void set_solution_vector(const NumericVector<double>& new_solution);


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


    //! Get the ID associated to a solution variable
    /*!
     *
     * \param solution_name the name of the solution variable
     * \return the ID, if \c solution_name exists, \c INVALID_ID otherwise
     */
    ID get_solution_id(const std::string& solution_name) const;


    //! Get the descriptor of a solution variable ID
    /*!
     *
     * \param id the id of the solution variable
     * \return the descriptor (which can have the ID \c INVALID_ID if
     *   \c id does not exist)
     */
    const SolutionDescriptor& get_solution_descriptor(ID id) const;


    //! Get the descriptor of a solution variable
    /*!
     *
     * \param solution_name the name of the solution variable
     * \return the descriptor (which can have the ID \c INVALID_ID if
     *   \c solution_name does not exist)
     */
    const SolutionDescriptor& get_solution_descriptor(const std::string& solution_name) const;


    //! Get solutions at specified points in an element
    /*!
     * \param elem the pointer to the element
     * \param values a map to hold the values, the key IDs specify the solutions
     *  to be returned
     * \param points the vector containing the points
     * \param local_coords set to true if the points are already in local coordinates
     * \return \c false if the model has not been solved or the element is invalid
     *
     * \note for efficiency it is preferable to pass local coordinates (i.e. coordinates
     *  in reference element)
     *
     * If a solution associated to a certain ID does not exist in the specified
     * element (e.g. if the element is not in the simulation domain), the corresponding
     * vector in the map has size zero. Generally the size of the vector depends on the
     * type of the solution.
     *
     *
     * If \c points is not specified, the values are returned on the ``natural''
     * locations of the solution variables.
     * \see SolutionDescriptor
     */
    bool get_solution(const Elem* elem, std::map<ID, std::vector<double> >& values,
        const std::vector<Point>& p = std::vector<Point>(),
        bool local_coords = false);

    //! Get a single solution
    /*!
     * \see get_solution(const Elem*, std::map<ID, std::vector<double> >&,
     *   const std::vector<Point>&, bool)
     */
    bool get_solution(const Elem* elem, ID id, std::vector<double>& values,
        const std::vector<Point>& p = std::vector<Point>(),
        bool local_coords = false);


    //! Get a single solution at a single point
    /*!
     * \see get_solution(const Elem*, std::map<ID, std::vector<double> >&,
     *   const std::vector<Point>&, bool)
     */
    bool get_solution(const Elem* elem, ID id, double& value,
        const Point& p,
        bool local_coords = false);


    //! Get solutions associated to an atom
    /*!
     *
     * \param atom a pointer to the atom
     * \param values a map to hold the values, the key IDs specify the solutions
     *  to be returned
     * \return \c false if the model has not been solved or the atom is invalid
     *
     * If a solution associated to a certain ID does not exist in the specified
     * element (e.g. if the element is not in the simulation domain), the corresponding
     * vector in the map has size zero. Generally the size of the vector depends on the
     * type of the solution (more specifically on the number of locations in the element
     * and on the number of solution components).
     * \see SolutionDescriptor
     */
    bool get_solution(const Atom* atom, std::map<ID, std::vector<double> >& values);


    //! Get solutions that are not associated to elements or atoms
    /*!
     *
     * \param values a map to hold the values, the key IDs specify the solutions
     *  to be returned
     * \return \c false if there is no solution
     *
     * If a solution associated to a certain ID does not exist, the corresponding
     * vector in the map has size zero. Generally the size of the vector depends on the
     * type of the solution (more specifically on the number of solution components).
     * \see SolutionDescriptor
     */
    bool get_solution(std::map<ID, std::vector<double> >& values);


    /*!
     * \copydoc build_nodal_results()
     *
     * Calls build_nodal_results()
     */
    void get_nodal_results(std::vector<double>& results,
        std::vector<std::string>& legend);


    /*!
     * \copydoc build_elemental_results()
     *
     * Calls build_elemental_results()
     */
    void get_elemental_results(std::vector<double>& results,
        std::vector<std::string>& legend);


    //! Get the type of this simulation
    /*!
     * The type is the identifying string which defines at creation time
     * which simulation to create. It's the same string one writes in the
     * input file.
     */
    const std::string& get_type(void) const;


    //! Set the scaling parameters
    void set_scaling(const Scaling& scaling);


    //! Get the scaling parameters
    Scaling& get_scaling(void);


    //! Get the mesh
    /*!
     * \return a constant reference to the simulation mesh
     */
    MeshBase& get_mesh(void) const;


    //! Get the mesh units (in meters)
    double get_mesh_units(void) const;


    //! Get the atomistic structure pointer
    /*!
     * Note that the given pointer may be NULL
     */
    AtomisticStructure* get_atomistic_structure(void) const;




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


    //! Tell the level of verbosity
    int verbose(void) const;


    //! Set the level of verbosity
    int& verbose(void);


    //! Check if \c region_id is include in this simulation
    bool includes_region(ID region_id) const;


    //! Create a bulk physical model to be used with this simulation
    PhysicalModel* new_bulk_model(const ModelOptions& options,
        const Material* material);


    //! Create a boundary model to be used with this simulation
    /*!
     * deprecated the use of BoundaryProperties is obsolete
     */
    BoundaryProperties* new_boundary_model(const ModelOptions& options);


    //! Create a boundary model to be used with this simulation
    PhysicalModel* new_boundary_model(const ModelOptions& options,
        const MaterialBoundary* boundary);


    //! Create an edge model
    PhysicalModel* new_edge_model(const ModelOptions& options,
        const EdgeObject* edge);


    //! Create a nodal model
    PhysicalModel* new_node_model(const ModelOptions& options,
        const NodeObject* node);


    //! Get the physical model for a certain region ID
    /*!
     * \return \c NULL if no model is present for region \c region_id
     * \deprecated get_bulk_model should be used instead
     */
    PhysicalModel* get_physical_model(ID region_id) const;


    //! Get the material for a given element
    Material* get_material(const Elem* elem) const;


    //! Get the physical model for a given element
    /*!
     * \return \c NULL if no model is present for \c elem
     */
    template <typename T>
    T* get_bulk_model(const Elem* elem) const;


    //! Get the physical model associated to an element side
    /*!
     * \return \c NULL if no model is present for the given side
     *
     * The model will have the bulk material set to the material
     * of \c elem
     */
    template <typename T>
    T* get_interface_model(const Elem* elem, int side) const;


    //! Get the physical model associated to an element edge
    /*!
     * \return \c NULL if no model is present for the given edge
     *
     * The model will have the bulk material set to the material
     * of \c elem
     */
    template <typename T>
    T* get_edge_model(const Elem* elem, int edge) const;


    //! Get the physical model associated to an element node
    /*!
     * \return \c NULL if no model is present for the given node
     *
     * The model will have the bulk material set to the material
     * of \c elem
     */
    template <typename T>
    T* get_node_model(const Elem* elem, int node) const;




    //! Get a reference to the set of all bulk physical models
    const std::set<PhysicalModel*>& get_physical_models(void) const;


    //! Get a reference to the set of all interface models
    const std::set<PhysicalModel*>& get_interface_models(void) const;


    //! Create an embracing region
    /*!
     * \param other_simulation the other simulation, the embracin region will
     * lie in that one.
     * \param options the options
     * \param need_mixing_coeff if \c true, the mixing coefficients
     * will be calculated, too
     */
    Embracing* create_embracing_region(SimulationInterface* other_simulation,
        const ModelOptions& options,
        bool need_mixing_coeff = false);


    //! Get the set of plotvariable IDs
    const IDSet& get_plotvariable_ids(void) const;



  protected:
 

    //! Empty constructor
    SimulationInterface(const ModelOptions& options);


    //! Get a reference to the equation system object
    EquationSystems& get_equation_systems(void) const;


    //! Clear systems
    void clear_systems(void);


    //! Get the unique name for the equation system
    const std::string& get_equation_system_name(void) const;


    //! Create a new equation system object
    /*!
     * \param type the type as string (linear, nonlinear, eigen)
     * The systems are stored in sequence of creation.
     */
    ID create_equation_system(const std::string& type);


    //! Get the equation system with number \c i
    /*!
     * \param i the index of the equation system
     * (the same as the sequence number at creation time)
     */
    template <typename T>
    T& get_equation_system(ID i = 0);


    //! Get the solver options
    ModelOptions& get_solver_options(void);


    //! Initialize the internal mesh pointer
    /*!
     * The default implementation just takes the mesh pointer from the device,
     * but this method may be reimplemented to do more exotic things.
     */
    virtual void setup_mesh(void);


    //! Initialize the internal atomistic structure pointer
    /*!
     * The default implementation just takes the structure pointer from the device,
     * getting the structure name from the input file ("atomistic_structure").
     * This method may be reimplemented to do more exotic things.
     * It may be good habit to reimplement it as empty for modules that cannot
     * do anything with atoms.
     */
    virtual void setup_atomistic_structure(void);


    //! Setup the available variables
    virtual void do_setup_solution_variables(void) {};


    //! Solve for equilibrium
    /*!
     * Can be implemented in derived classes to solve for some
     * equilibrium state.
     * It has to be called explicitly.
     *
     * Use it e.g. to solve the Poisson equation for equilibrium, to set
     * some reasonable starting values etc.
     */
    virtual void do_equilibrium(void) {};


    //! Do the initialization
    /*!
     * Has to be implemented by derived classes.
     *
     * This method should initialize everything that is needed to do a
     * simulation.
     */
    virtual void do_init(void) = 0;


    //! Do some reinitialization before any solve
    /*!
     * May be reimplemented by derived classes if necessary.
     */
    virtual void do_reinit(void) {};


    //! Do the solve
    /*!
     * Has to be implemented by derived classes.
     *
     * This method does the actual simulation.
     */
    virtual void do_solve(void) = 0;


    //! Increment the solve sequence number
    void increment_solve_sequence_number(void);


    //! Declare this module to be a task or not
    void is_task(bool task);


    //! Adds a solution name to the plot list
    /*!
     * This needs to be done \em before calling
     * declare_solution()
     */
    void add_plot_variable(const std::string& name);


    //! Adds a solution id to the plot variable list
    void add_plot_variable(ID id);


    //! Checks if a solution variable should be plotted
    bool plot_solution(const std::string& name) const;


    //! Checks if a solution variable with given ID should be plotted
    /*!
     * This method can be used only \em after init()
     */
    bool plot_solution(ID id) const;


    //! Do a plot
    /*!
     * This method creates files with the results of the simulation.
     * In almost all cases the default implementation should be ok, but
     * in some it is not, especially in "compound" simulations as
     * sweeps, selfconsitent solvers etc.
     *
     * It calls plot_meshdata(), plot_atomisticdata() and plot_globaldata()
     */
    virtual void do_plot(void);
    virtual void do_plot_old(void);


    //! Plot mesh associated data
    /*!
     * This method plots data associated with the simulation mesh. Usually,
     * the default implementation should be ok.
     */
    virtual void plot_meshdata(void);


    //! Plot data associated with atoms
    /*!
     * This method plots data associated with the atoms mesh. Usually,
     * the default implementation should be ok.
     *
     * TODO this has to be implemented
     */
    virtual void plot_atomisticdata(void);


    //! Plot global data
    /*!
     * This method writes global data without association with a grid
     * (e.g. contact currents, energy levels etc.)
     *
     * The default implementation may not make much sense in certain
     * cases and should then be reimplemented.
     */
    virtual void plot_globaldata(void);


    //! Print simulation info
    virtual void do_print_info(void);


    //! Remember the current solution
    /*!
     * The default action is to clone the solution vector of the
     * equations system. If \c id is invalid, a new ID will be generated
     * and returned, if not, the remembered solution \c id will be overwritten.
     * If this method is reimplemented in a derived class, it should behave in
     * the same way!
     */
    virtual ID do_remember_current_solution(ID id = 0);


    //! Tells the base class if the model has a solution vector
    /*!
     * \param flag = true: the model has a solution vector
     *
     * Call this method in the constructor if your model does \em not
     * have a solution vector.
     */
    void has_solution_vector(bool flag);


    //! Get a pointer to the solution vector
    virtual NumericVector<double>& do_get_solution_vector(void);


    //! Set a new solution vector
    /*!
     * The default implementation just makes
     * \c solution_vector() = \c new_solution
     *
     * No checks will be done on the vector size!
     * \param new_solution the new solution to set
     */
    virtual void do_set_solution_vector(
        const NumericVector<double>& new_solution);


    //! Set to the remembered solution number \c id
    virtual void do_set_to_remembered_solution(ID id);


    //! Delete a remembered solution
    virtual void do_delete_remembered_solution(ID id);


    //! Build the maximum norm of the solution difference
    /*!
     * This method returns the maximum norm
     * \f$\Vert x - x_i\Vert_\infty\f$
     * of the difference between the current solution \f$x\f$ and some
     * remembered solution \f$x_i\f$
     *
     * \param id the id of the remembered solution.
     *
     * If the remembered solution doesn't exist, return value is zero.
     *
     * The default action should be ok in most cases. But it could
     * be useful to reimplement this method when a simulation
     * uses some scaling
     */
    virtual double do_maximum_norm_of_difference(ID id);


    //! Build the l2 norm of the solution difference
    /*!
     * This method returns the l2 norm \f$\Vert x - x_i\Vert_2\f$
     * of the difference between the current solution \f$x\f$ and some
     * remembered solution \f$x_i\f$
     *
     * \param id the id of the remembered solution.
     *
     * If the remembered solution doesn't exist, return value is zero.
     *
     * The default action should be ok in most cases. But it could
     * be useful to reimplement this method when a simulation
     * uses some scaling
     */
    virtual double do_l2_norm_of_difference(ID id);



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


    //! Get the output data formats
    /*!
     * At least one format will always be returned.
     */
    void get_output_format(std::vector<std::string>& formats) const;


    //! Get the output directory
    std::string get_output_directory(void) const;


    //! Get the output filename prefix
    std::string get_output_filename_prefix(void) const;


    //! Get the whole output filename (except extension)
    /*!
     * \return get_output_filename_prefix() +
     *    TiberCad::get_filename_suffix()
     */
    std::string get_output_filename(void) const;


    //! Check if binary data should be written
    bool binary_output(void) const;


    //! Get the ID for a given variable name
    /*!
     * An ID of \c INVALID_ID represents an unknown variable.
     *
     * \param the name of the variable as a string
     * \return the ID of the variable as a numerical value
     */
    virtual ID convert_variable_name_to_id(const std::string& variable_name) const;


    //! Register a solution
    /*!
     *
     * \param name the name of the solution
     * \param id the ID to be assigned to the solution
     * \param type the type of the solution (see SolutionDescriptor)
     * \param location the location of the quantity inside an element
     *        (see SolutionDescriptor)
     * \param units the physical units, if applicable (recommended)
     * \param n_comp the number of components (only needed for generic n-tuples)
     *
     * It is convenient to use this method via the declare_solution macro. This
     * macro accepts as the first argument an enum value (or the name of a static
     * member variable) whose name as a string will be used as the solution
     * variable name. For \c type and \c location one can use directly the enum values (instead
     * of the fully qualified name).
     * Example:
     * \code
     * enum Solutions {
     *   Field,
     *   Potential
     * };
     *
     * ...
     *
     * declare_solution(Field, VECTOR, CELL, "V/m");
     * declare_solution(Potential, REAL, NODAL, "V");
     * \endcode
     *
     * \note You can call declare_solution() more than once for the same id, but this
     * will overwrite the solution descriptor inserted before!
     */
    void declare_solution_ext(const std::string& name, ID id,
        SolutionDescriptor::Type type, SolutionDescriptor::Location location,
        const std::string& units = "", unsigned int n_comp = 0);

    /*!
     * \brief This is a shorter version of declare_solution_ext()
     *
     * \li \c name the name of the solution
     * \li \c type the type of the solution (see SolutionDescriptor)
     * \li \c location the location of the quantity inside an element
     *        (see SolutionDescriptor)
     * \li \c units the physical units, if applicable (optional, recommended)
     * \li \c n_comp the number of components (only needed for generic n-tuples)
     *
     * The name should correspond to an enum value (or the name of a static
     * member variable).
     * For \c type and \c location one can use directly the enum values (instead
     * of the fully qualified name)
     *
     * \protected \memberof SimulationInterface
     */
#define declare_solution(name, type, location, ...) \
    declare_solution_ext(#name, name, SolutionDescriptor::type, \
        SolutionDescriptor::location, ## __VA_ARGS__)

    //! Add an alias for a solution variable
    void add_alias(const std::string& alias, ID id);



    //! Get solutions at specified points in an element
    /*!
     * \param elem the pointer to the element
     * \param values a map to hold the values, the key IDs specify the solutions
     *  to be returned
     * \param p the vector containing the points
     *
     * \note The points in \c p are given in \em local coordinates (i.e. on the
     * reference element)
     *
     * \note Cell based solutions have to be returned only once, not for every
     * requested point.
     *
     * \pre \c elem is an active element of this simulation
     * \pre all ids refer to solutions located on an element
     */
    virtual void get_solution_secure(const Elem* elem,
        std::map<ID, std::vector<double> >& values,
        const std::vector<Point>& p);


    //! Get solutions on an atom
    /*!
     * \param atom the pointer to the atom
     * \param values a map to hold the values, the key IDs specify the solutions
     *  to be returned
     *
     * \pre \c atom is an atom of this simulation
     * \pre all ids refer to solutions associated to an atom
     */
    virtual void get_solution_secure(const Atom* atom,
        std::map<ID, std::vector<double> >& values);


    //! Get solutions not located on the mesh or the atoms
    /*!
     * \param values a map to hold the values, the key IDs specify the solutions
     *  to be returned
     *
     *  \pre all ids refer to mesh independent solutions
     */
    virtual void get_solution_secure(std::map<ID, std::vector<double> >& values);


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


    //! Create a physical model that can be used by this type of simulation
    /*!
     * The default behaviour defined in the base class is to return the
     * NULL pointer, because there could be simulations that don't need
     * a physical model. If a derived class reimplements this method (which
     * will normally be the case) it should notify about errors by throwing
     * a ModelErrorException.
     *
     * \param options the options
     * \param mat the material this model is associated with
     */
    virtual PhysicalModel* create_bulk_model(const ModelOptions& options,
        const Material* mat) const;

    virtual PhysicalModel*
      create_physical_model(const ModelOptions& options,
                            const Material* mat) const
      throw (ModelErrorException);


    //! Create a boundary model that can be used by this type of simulation
    /*!
     * The default behaviour defined in the base class is to return the
     * NULL pointer, because there could be simulations that don't need
     * a boundary model. If a derived class reimplements this method (which
     * will normally be the case) it should notify about errors by throwing
     * a ModelErrorException.
     *
     * \param options the options
     * \param boundary the boundary object this model is associated with
     */
    virtual PhysicalModel* create_boundary_model(const ModelOptions& options,
        const MaterialBoundary* boundary) const;

    virtual BoundaryProperties*
      create_boundary_model(const ModelOptions& options) const
      throw (ModelErrorException);


    //! Create an edge model that can be used by this type of simulation
    /*!
     * The default behaviour defined in the base class is to return the
     * NULL pointer, because there could be simulations that don't need
     * an edge model. If a derived class reimplements this method (which
     * will normally be the case) it should notify about errors by throwing
     * a ModelErrorException.
     *
     * \param options the options
     * \param edge the edge object this model is associated with
     */
    virtual PhysicalModel* create_edge_model(const ModelOptions& options,
        const EdgeObject* edge) const;


    //! Create a nodal model that can be used by this type of simulation
    /*!
     * The default behaviour defined in the base class is to return the
     * NULL pointer, because there could be simulations that don't need
     * an nodal model. If a derived class reimplements this method (which
     * will normally be the case) it should notify about errors by throwing
     * a ModelErrorException.
     *
     * \param options the options
     * \param node the node object this model is associated with
     */
    virtual PhysicalModel* create_node_model(const ModelOptions& options,
        const NodeObject* node) const;


    //! Set or unset the \c equilibrium_done flag
    void equilibrium_done(bool flag);


    //! Save the current state
    /*!
     * This method may be overridden. The default implementation
     * writes the vector obtained from get_solution_vector()
     * (if available) into a file.
     */
    virtual void do_save_data(std::ostream& os);


    //! Reload the current state
    /*!
     * This method may be overridden. The default implementation
     * reads data into the vector obtained from get_solution_vector()
     * (if available).
     */
    virtual void do_load_data(std::istream& is);



  private:

    //! A typedef for the embracing region map
    typedef std::map<SimulationInterface*, Embracing*> EmbracingMap;


    //! The type of the solution descriptor map
    typedef std::map<ID, SolutionDescriptor> SolutionDescrMap;


    //! Some common solution variables
    /*!
     * We assign the IDs below INVALID_ID so a clash with the ones
     * defined in modules should be impossible.
     */
    enum Solution
    {
      RegionIDs = (INVALID_ID - 1)        /*!< The region IDs */
    };

    //! The environment for this simulation
    SimulationEnvironment* _environment;


    //! A flag indicating if the simulator is initialized
    bool _is_initialized;


    //! The solve sequence number
    unsigned int _solve_sequence_nr;


    //! \c true if this module is a task module
    bool _is_task;


    //! A flag indicating that equilibrium has been done
    bool _equilibrium_is_solved;


    //! Do we have a solution vector or not
    bool _has_solution_vector;



    //! The ID of this simulation
    /*!
     * The ID is unique for every simulator and is assigned automatically
     * at instantiation.
     */
    ID _id;


    //! The identifying string for the type of this simulation
    std::string _type;


    //! The unique name for the equation system
    std::string _eq_system_name;


    //! The scaling parameters
    Scaling _scaling;


    //! A map with all embracing regions
    EmbracingMap _embracings;


    //! A map with the descriptors of all known solutions
    SolutionDescrMap _solution_descriptors;


    //! A map mapping solution names to IDs
    std::map<const std::string, ID> _solution_ids;


    //! A map with remembered solutions
    std::map<ID, NumericVector<double>*> _remembered_solutions;


    //! A set with all physical models of this simulation
    std::set<PhysicalModel*> _physical_models;


    //! A set with all boundary models of this simulation
    std::set<PhysicalModel*> _boundary_models;


    //! The level of verbosity
    int _verbosity;


    //! The set of all variables to be plotted
    std::set<std::string> _plotvariables;


    //! The set of all variable IDs to be plotted
    IDSet _plotvariable_ids;


    //! The map containing all simulations with their ID
    static SimulationMap _simulation_map;


    //! An invalid solution descriptor
    /*!
     * Some accessor methods need to return a reference to
     * an invalid solution descriptor.
     */
    static SolutionDescriptor _invalid_descr;


    //! The equation systems
    std::vector<TiberEqSystem*> _systems;


    //! A pointer to the simulation mesh
    MeshBase* _mesh;


    //! A pointer to the principal atomistic structure
    AtomisticStructure* _atomistic_structure;


    //! create a unique name for the equation system
    void create_equation_system_name(void) TBDLLOCAL;


    //! Set the simulation type (= identifier)
    /*!
     * The identifier is used at creation time to know which type of
     * simulation to create.
     */
    void set_type(const std::string& type) TBDLLOCAL;


    //! Print all registered solution variables
    void print_known_solution_variables(void) const TBDLLOCAL;


    //! Do not allow copy constructor
    SimulationInterface(const SimulationInterface&) TBDLLOCAL;

    //! Do not allow assignement operator
    SimulationInterface& operator=(const SimulationInterface&) TBDLLOCAL;


    //! \see get_bulk_model()
    PhysicalModel* _get_bulk_model(const Elem* elem) const;


    //! \see get_interface_model()
    PhysicalModel* _get_interface_model(const Elem* elem, int side) const;


    //! \see get_edge_model()
    PhysicalModel* _get_edge_model(const Elem* elem, int edge) const;


    //! \see get_node_model()
    PhysicalModel* _get_node_model(const Elem* elem, int node) const;
};


//
// inline members
//




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


template <typename T>
inline
T&
SimulationInterface::get_equation_system(ID i)
{
  if (i >= _systems.size())
    throw RuntimeException("Trying to access inexistent system.");

#ifdef DEBUG
  return dynamic_cast<T&>(*_systems[i]);
#else
  return static_cast<T&>(*_systems[i]);
#endif
}





inline
bool
SimulationInterface::has_environment(void) const
{
  return (_environment != 0);
}


inline
SimulationEnvironment&
SimulationInterface::get_environment(void) const
{
  if (_environment == NULL)
    throw RuntimeException("Simulation \'" + get_name() +
        "\' has no simulation environment associated");

  return *_environment;
}


inline
MeshBase&
SimulationInterface::get_mesh(void) const
{
  if (_environment == NULL)
    throw RuntimeException("Simulation \'" + get_name() +
        "\' has no mesh associated");

  return *_mesh;
}



inline
AtomisticStructure*
SimulationInterface::get_atomistic_structure(void) const
{
  return _atomistic_structure;
}


inline
const std::set<PhysicalModel*>&
SimulationInterface::get_physical_models(void) const
{
  return _physical_models;
}



inline
const std::set<PhysicalModel*>&
SimulationInterface::get_interface_models(void) const
{
  return _boundary_models;
}



inline
ID
SimulationInterface::get_id(void) const
{
  return _id;
}




inline
const std::string&
SimulationInterface::get_equation_system_name(void) const
{
  return _eq_system_name;
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
  return (get_solve_sequence_number() > 0);
}




inline
bool
SimulationInterface::is_task(void) const
{
  return _is_task;
}




inline
void
SimulationInterface::is_task(bool task)
{
  _is_task = task;
}


inline
const IDSet&
SimulationInterface::get_plotvariable_ids(void) const
{
  return _plotvariable_ids;
}


inline
bool
SimulationInterface::plot_solution(const std::string& name) const
{
  return _plotvariables.count(name);
}


inline
bool
SimulationInterface::plot_solution(ID id) const
{
  return _plotvariable_ids.count(id);
}



inline
unsigned int
SimulationInterface::get_solve_sequence_number(void) const
{
  return _solve_sequence_nr;
}


inline
bool
SimulationInterface::equilibrium_done(void) const
{
  return _equilibrium_is_solved;
}



inline
void
SimulationInterface::equilibrium_done(bool flag)
{
  _equilibrium_is_solved = flag;
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
SimulationInterface::build_nodal_results(
    const std::set<std::string>&,
    std::vector<double>&, std::vector<std::string>&)
{
}



inline
void
SimulationInterface::build_elemental_results(
    const std::set<std::string>&,
    std::vector<double>&, std::vector<std::string>&)
{
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
SimulationInterface::get_l2_norm_of_difference(ID id)
{
  return do_l2_norm_of_difference(id);
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
bool
SimulationInterface::get_solution(const Elem* elem, ID id,
    double& value, const Point& p,
    bool local_coords)
{
  std::vector<Point> point(1, p);
  std::vector<double> val(1, 0.0);

  bool ok = get_solution(elem, id, val, point);
  if (ok)
    value = val[0];

  return ok;
}

/*
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
*/

inline
NumericVector<double>&
SimulationInterface::get_solution_vector(void)
{
  assert(has_solution_vector());

  return do_get_solution_vector();
}


inline
void
SimulationInterface::set_solution_vector(const NumericVector<double>& new_solution)
{
  if (has_solution_vector())
    do_set_solution_vector(new_solution);
}


inline
bool
SimulationInterface::has_solution_vector(void)
{
  return _has_solution_vector;
}


inline
void
SimulationInterface::has_solution_vector(bool flag)
{
  _has_solution_vector = flag;
}



inline
int
SimulationInterface::verbose(void) const
{
  return _verbosity;
}



inline
int&
SimulationInterface::verbose(void)
{
  return _verbosity;
}


template <typename T>
inline
T*
SimulationInterface::get_bulk_model(const Elem* elem) const
{
  return dynamic_cast<T*>(_get_bulk_model(elem));
}


template <typename T>
inline
T*
SimulationInterface::get_interface_model(const Elem* elem, int side) const
{
  return dynamic_cast<T*>(_get_interface_model(elem, side));
}


template <typename T>
inline
T*
SimulationInterface::get_edge_model(const Elem* elem, int edge) const
{
  return dynamic_cast<T*>(_get_edge_model(elem, edge));
}



template <typename T>
inline
T*
SimulationInterface::get_node_model(const Elem* elem, int node) const
{
  return dynamic_cast<T*>(_get_node_model(elem, node));
}


inline
SimulationInterface::SimulationIterator
SimulationInterface::simulations_begin(void)
{
  return SimulationIterator(_simulation_map.begin());
}


inline
SimulationInterface::SimulationIterator
SimulationInterface::simulations_end(void)
{
  return SimulationIterator();
}

inline
void
SimulationInterface::increment_solve_sequence_number(void)
{
  ++_solve_sequence_nr;
}

#endif // _SIMULATIONINTERFACE_H_
