// $Id$

#ifndef _SELFCONSISTENTSOLVER_H_
#define _SELFCONSISTENTSOLVER_H_

#include "SimulationInterface.h"


class Multiscale;
class XMonitor;

//! Interface for self-consistent calculations
class TBDLLOCAL SelfconsistentSolver : public SimulationInterface
{

  public:

    //! The destructor
    virtual ~SelfconsistentSolver();


  protected:

    //! An iterator for the simulations
    typedef std::vector<SimulationInterface*>::iterator Iterator;

    //! The empty Constructor
    SelfconsistentSolver(const ModelOptions& options);

    
    /*! \copydoc SimulationInterface::do_init() */
    virtual void do_init(void);

    
    /*! \copydoc SimulationInterface::do_equilibrium() */
    virtual void do_equilibrium(void);

    
    /*! \copydoc SimulationInterface::do_solve() */
    virtual void do_solve(void) = 0;

    
    /*! \copydoc SimulationInterface::do_plot() */
    virtual void do_plot(void);

    
    /*! \copydoc SimulationInterface::parse_options()
     * If this method is reimplemented in a derived class
     * it has to called explicitly by \c SelfconsistentSolver::parse_options()
     * */
    virtual void parse_options(void);


    /*! \copydoc SimulationInterface::do_remember_current_solution() */
    virtual ID do_remember_current_solution(ID id = 0);


    /*! \copydoc SimulationInterface::do_set_to_remembered_solution() */
    virtual void do_set_to_remembered_solution(ID id);


    /*! \copydoc SimulationInterface::do_delete_remembered_solution() */
    virtual void do_delete_remembered_solution(ID id);
    

    //! Get global quantities from all sub-simulations
    virtual void get_solution_secure(std::map<ID, std::vector<double> >& values);


    //! Returns the solution vector of the last simulation
    virtual NumericVector<double>& do_get_solution_vector(void);


    //! Returns the solution vector of the last simulation
    virtual void do_set_solution_vector(const NumericVector<double>& new_solution);


    //! Abused here for the inheritance of variables
    virtual void do_print_info(void);

    
    //! Initialize the solver
    /*!
     * Can be used to initialize the solutions.
     *
     * The default implementation just solves once all simulations.
     */
    virtual void initialize(void);


    //! Solve all simulations
    void solve_simulations(void);


    //! Get the maximum number of iterations
    unsigned int get_maximum_iterations(void) const;


    //! Get the relative tolerance
    double get_relative_tolerance(void) const;


    //! Get the absolute tolerance
    double get_absolute_tolerance(void) const;


    //! Get monitor
    bool get_monitor(void) const;


    //! Get the pointer to the X monitor
    /*!
     * Use this with caution, especially don't forget to check for NULL pointer
     */
    XMonitor* get_xmonitor(void);


    //! Opens the X monitor
    /*!
     * Call this in do_solve() at the beginning
     */
    void open_xmonitor(void);


    //! Closes the X monitor
    /*!
     * Call this at the end of do_solve()
     */
    void close_xmonitor(void);


    //! Add a point to the X monitor
    /*!
     * \param iteration the iteration number
     * \param err the error
     * \param logarithm if \c true, plot \c log10(error)
     */
    void draw_point(double iteration, double error, bool logarithm = true);



    //! Get an iterator for the first simulation
    Iterator simulations_begin(void);


    //! Get the past-the-end iterator for the simulations
    Iterator simulations_end(void);


    //! Get the number of simulations
    int get_number_of_simulations(void) const;


    //! Get the i-th simulation
    /*!
     * Returns \c NULL if index \c i >= # of simulations
     */
    SimulationInterface* simulation(unsigned int i);


    //! Get the last simulation, that provides also the solution vector
    SimulationInterface* get_last_simulation(void);


  private:

    //! The simulations which we solve self-consistently
    std::vector<SimulationInterface*> _simulations;

    //! The maximum number of iterations
    unsigned int _max_it;

    //! The relative tolerance
    double _rel_tol;

    //! The absolute tolerance
    double _abs_tol;

    //! The ids of the remembered solutions
    std::map<ID, std::vector<ID> > _remembered_sol_ids;

    //! Screen output of the convergence process 
    bool _monitor;
    
    //! The X monitor
    XMonitor* _xmonitor;


    //! Multiscale definitions
    Multiscale* _multiscale;

};


//
// inline methods
//

inline
SelfconsistentSolver::SelfconsistentSolver(const ModelOptions& options)
  : SimulationInterface(options),
    _max_it(20),
    _rel_tol(1e-3),
    _abs_tol(1e-3),
    _xmonitor(NULL)
{
  is_task(true);
}



inline
unsigned int
SelfconsistentSolver::get_maximum_iterations(void) const
{
  return _max_it;
}

inline
double
SelfconsistentSolver::get_relative_tolerance(void) const
{
  return _rel_tol;
}

inline
double
SelfconsistentSolver::get_absolute_tolerance(void) const
{
  return _abs_tol;
}


inline 
bool
SelfconsistentSolver::get_monitor(void) const
{
  return _monitor;
}


inline 
XMonitor*
SelfconsistentSolver::get_xmonitor(void)
{
  return _xmonitor;
}



inline
SelfconsistentSolver::Iterator
SelfconsistentSolver::simulations_begin(void)
{
  return _simulations.begin();
}

inline
SelfconsistentSolver::Iterator
SelfconsistentSolver::simulations_end(void)
{
  return _simulations.end();
}


inline
int
SelfconsistentSolver::get_number_of_simulations(void) const
{
  return _simulations.size();
}


inline
SimulationInterface*
SelfconsistentSolver::simulation(unsigned int i)
{
  if (i < _simulations.size())
    return _simulations[i];
  else
    return NULL;
}


inline
SimulationInterface*
SelfconsistentSolver::get_last_simulation(void)
{
  return _simulations[_simulations.size() - 1];
}




#endif // _SELFCONSISTENTSOLVER_H_
