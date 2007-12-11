// $Id$

#ifndef _SELFCONSISTENTSOLVER_H_
#define _SELFCONSISTENTSOLVER_H_

#include "SimulationInterface.h"


//! Interface for self-consistent calculations
class SelfconsistentSolver : public SimulationInterface
{

  public:

    //! The destructor
    virtual ~SelfconsistentSolver(void);


  protected:

    //! An iterator for the simulations
    typedef std::vector<SimulationInterface*>::iterator SimulationIterator;

    //! The empty Constructor
    SelfconsistentSolver(void);

    
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
    

    /*! \copydoc SimulationInterface::build_integrated_quantities() */
    virtual void build_integrated_quantities(
        const std::set<std::string>& names,
        std::vector<double>& values);


    /*! \copydoc SimulationInterface::build_integrated_quantities() */
    virtual void build_integrated_quantities_description(
        const std::set<std::string>& names,
        std::vector<std::string>& legend,
        std::vector<std::string>& description);


    //! Returns the solution vector of the last simulation
    virtual NumericVector<double>& do_get_solution_vector(void);


    //! Returns the solution vector of the last simulation
    virtual void do_set_solution_vector(const NumericVector<double>& new_solution);

    
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


    //! Get xmonitor
    bool get_xmonitor(void) const;

    //! Get an iterator for the first simulation
    SimulationIterator simulations_begin(void);


    //! Get the past-the-end iterator for the simulations
    SimulationIterator simulations_end(void);


    //! Get the number of simulations
    int get_number_of_simulations(void) const;


    //! Get the i-th simulation
    /*!
     * Returns \c NULL if index \c i >= # of simulations
     */
    SimulationInterface* simulation(unsigned int i);



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

    //! screen output of the convergence process 
    bool _monitor;
    
    //! graphical output of the convergence process 
    bool _xmonitor;

};


//
// inline methods
//

inline
SelfconsistentSolver::SelfconsistentSolver(void)
  : _max_it(5),
    _rel_tol(1e-3),
    _abs_tol(1e-3)
{
}

inline
SelfconsistentSolver::~SelfconsistentSolver(void)
{
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
bool
SelfconsistentSolver::get_xmonitor(void) const
{
  return _xmonitor;
}


inline
SelfconsistentSolver::SimulationIterator
SelfconsistentSolver::simulations_begin(void)
{
  return _simulations.begin();
}

inline
SelfconsistentSolver::SimulationIterator
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




#endif // _SELFCONSISTENTSOLVER_H_
