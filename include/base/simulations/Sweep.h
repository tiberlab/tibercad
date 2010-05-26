// $Id$

#ifndef _SWEEP_H_
#define _SWEEP_H_

#include "SimulationInterface.h"

#include <set>
#include <vector>



//! Make a sweep of a sweepable variable
/*!
 * This class takes a variable (cf. Variable class)
 * and performs a sweep with it (e.g. contact voltage). For each
 * sweep step it solves the given simulation(s) according to
 * the following flow chart:
 * \image html Sweep_flowchart.jpg
 * \image latex Sweep_flowchart.eps
 *    "Flow chart for a parameter sweep" width=10cm
 */
class Sweep : public SimulationInterface
{

  public:
    
    //! Destructor
    virtual ~Sweep(void);

    //! Create a Sweep object
    static Sweep* create(const ModelOptions& options);


    //! Get the list of inner simulations
    const std::vector<SimulationInterface*>
      get_simulations(void) const;



  protected:

    //! Constructor
    Sweep(const ModelOptions& options);

    /*! \copydoc SimulationInterface::do_init() */
    virtual void do_init(void);


    /*! \copydoc SimulationInterface::do_equilibrium() */
    virtual void do_equilibrium(void);

    
    /*! \copydoc SimulationInterface::do_solve()
     *
     * solve() of the associated SimulationInterface objects will
     * be called repeatedly for all sweep values
     */
    virtual void do_solve(void);
   
    //! Plot results. Empty implementation!
    /*!
     * We do not plot anything here because we plot after each sweep step
     */
    virtual void do_plot(void);


    /*! \copydoc SimulationInterface::do_remember_current_solution() */
    virtual ID do_remember_current_solution(ID id = 0);


    /*! \copydoc SimulationInterface::do_set_to_remembered_solution() */
    virtual void do_set_to_remembered_solution(ID id);


    /*! \copydoc SimulationInterface::do_delete_remembered_solution() */
    virtual void do_delete_remembered_solution(ID id);
  

    
    /*! \copydoc SimulationInterface::parse_options() */
    virtual void parse_options(void);
    

    //! Get global quantities from all sub-simulations
    virtual void get_solution_secure(std::map<ID, std::vector<double> >& values);


    //! Get the innermost simulations
    /*!
     * Expands all sweeps, not including the simulations following
     * after a nested sweep.
     */
    void get_inner_simulations(
        std::vector<SimulationInterface*>& sims) const;


  private:

    //! The simulations for which wew do the sweep
    std::vector<SimulationInterface*> _simulations;

    //! A pointer to the sweepable variable
    std::string _variable;

    //! A vector containing all the simulation values
    std::vector<double> _values;

    //! The current simulation goal
    double _goal;

    //! The minimum step size
    double _min_step;

    //! The maximum step size
    double _max_step;


    //! Write results to file after every step if true
    bool _plot_data;


    //! The current filename suffix
    std::string _suffix;


    //! The ids of the remembered solutions
    std::map<ID, std::vector<ID> > _remembered_sol_ids;

    //! The remembered sweep value
    double _remembered_sweep_value;


    //! Remember the current solution
    ID remember_solution(void);


    //! Do the sweep
    void do_sweep(std::vector<double>& values,
        std::vector<std::ofstream*>& plotfiles,
        std::vector<std::map<double, std::vector<double> > >& sweep_data);


    //! Write the global data to file
    void write_global_data(SimulationInterface& simulation, std::ofstream*& file);


};


//
// inline methods
//


inline
Sweep::Sweep(const ModelOptions& options)
  : SimulationInterface(options),
    _variable(""),
    _min_step(1e-3),
    _max_step(10),
    _plot_data(false)
{
  has_solution_vector(false);
}


inline
Sweep*
Sweep::create(const ModelOptions& options)
{
  return new Sweep(options);
}


inline
const std::vector<SimulationInterface*>
Sweep::get_simulations(void) const
{
  return _simulations;
}


#endif // _SWEEP_H_
