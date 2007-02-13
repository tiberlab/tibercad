// $Id$

#ifndef _SWEEP_h_
#define _SWEEP_h_

#include "SimulationInterface.h"

#include <set>


class Sweepable;

//! Make a sweep of a sweepable variable
/*!
 * This class takes a sweepable (cf. Sweepable class) boundary model
 * and makes a sweep over that variable (e.g. contact voltage). For each
 * sweep step it solves the given simulation according to the following
 * flow chart:
 * \image html Sweep_flowchart.jpg
 * \image latex Sweep_flowchart.eps "Flow chart for a parameter sweep" width=10cm
 */
class Sweep : public SimulationInterface
{

  public:

    //! Constructor
    Sweep(void);
    
    //! Destructor
    virtual ~Sweep(void);

    //! Create a Sweep object
    static Sweep* create(void);


  protected:

    /*! \copydoc SimulationInterface::do_init() */
    virtual void do_init(void);

    /*! \copydoc SimulationInterface::do_solve()
     * solve() of the associated SimulationInterface object will
     * be called repeatedly for all sweep values
     */
    virtual void do_solve(void);

    /*! \copydoc SimulationInterface::parse_options() */
    virtual void parse_options(void);
    

  private:

    //! The simulation for which wew do the sweep
    SimulationInterface* _simulation;

    //! A pointer to the sweepable variable
    Sweepable* _variable;

    //! A vector containing all the simulation values
    std::vector<double> _values;

    //! The last simulation value
    double _last;

    //! The minimum step size
    double _min_step;

    //! The maximum step size
    double _max_step;

    //! A pointer to a vector containing always the last solution
    NumericVector<Real>* _last_solution;

    //! Write results to file after every step if true
    bool _do_output;


    //! The dependent variables we want to plot
    std::set<std::string> _plotvariables;


};


//
// inline methods
//


inline
Sweep::Sweep(void)
  : _simulation(0),
    _variable(0),
    _min_step(1e-3),
    _max_step(1),
    _last_solution(0),
    _do_output(true)
{
}


inline
Sweep::~Sweep(void)
{
}


inline
Sweep*
Sweep::create(void)
{
  return new Sweep();
}


#endif // _SWEEP_h_
