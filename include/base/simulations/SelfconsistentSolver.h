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

    //! Create a Sweep object
    static SelfconsistentSolver* create(void);


  protected:

    //! The empty Constructor
    SelfconsistentSolver(void);

    
    /*! \copydoc SimulationInterface::do_init() */
    virtual void do_init(void);

    
    /*! \copydoc SimulationInterface::do_equilibrium() */
    virtual void do_equilibrium(void);

    
    /*! \copydoc SimulationInterface::do_solve() */
    virtual void do_solve(void);

    
    /*! \copydoc SimulationInterface::do_plot() */
    virtual void do_plot(void);

    
    /*! \copydoc SimulationInterface::parse_options() */
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



  private:

    //! The simulations which we solve self-consistently
    std::vector<SimulationInterface*> _simulations;

    //! The maximum number of iterations
    unsigned int _max_it;

    //! The relative tolerance
    double _rel_tol;

    //! The absolute tolerance
    double _abs_tol;

    //! The relaxation factor to be used
    double _relax;

    //! The ids of the remembered solutions
    std::map<ID, std::vector<ID> > _remembered_sol_ids;

};


//
// inline methods
//

inline
SelfconsistentSolver::SelfconsistentSolver(void)
  : _max_it(5),
    _rel_tol(1e-3),
    _abs_tol(1e-3),
    _relax(0.5)
{
}

inline
SelfconsistentSolver::~SelfconsistentSolver(void)
{
}

inline
SelfconsistentSolver*
SelfconsistentSolver::create(void)
{
  return new SelfconsistentSolver();
}



#endif // _SELFCONSISTENTSOLVER_H_
