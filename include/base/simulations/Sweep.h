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
 * sweep step it solves the given simulation(s) according to the following
 * flow chart:
 * \image html Sweep_flowchart.jpg
 * \image latex Sweep_flowchart.eps "Flow chart for a parameter sweep" width=10cm
 */
class Sweep : public SimulationInterface
{

  public:
    
    //! Destructor
    virtual ~Sweep(void);

    //! Create a Sweep object
    static Sweep* create(void);


  protected:

    //! Constructor
    Sweep(void);

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
    

  private:

    //! The simulations for which wew do the sweep
    std::vector<SimulationInterface*> _simulations;

    //! A pointer to the sweepable variable
    std::string _variable;

    //! A vector containing all the simulation values
    std::vector<double> _values;

    //! The last simulation value
    double _last;

    //! The minimum step size
    double _min_step;

    //! The maximum step size
    double _max_step;


    //! Write results to file after every step if true
    bool _do_output;


    //! The dependent variables we want to plot
    std::set<std::string> _plotvariables;


    //! The ids of the remembered solutions
    std::map<ID, std::vector<ID> > _remembered_sol_ids;


    //! Remember the current solution
    ID remember_solution(void);


    //! Do the sweep
    void do_sweep(std::vector<double>& values,
        std::vector<std::ofstream*>& plotfiles,
        std::vector<std::map<double, std::vector<double> > >& sweep_data);


    //! Plot the data
    void plot_data(std::vector<std::ofstream*>& plotfiles,
        std::vector<std::map<double, std::vector<double> > >& sweep_data);

};


//
// inline methods
//


inline
Sweep::Sweep(void)
  : _variable(""),
    _min_step(1e-3),
    _max_step(1),
    _do_output(true)
{
}


inline
Sweep*
Sweep::create(void)
{
  return new Sweep();
}


#endif // _SWEEP_H_
