// $Id$

#ifndef _RAMP_H_
#define _RAMP_H_

#include "tiber_dll.h"
#include "TypeDefs.h"

#include <vector>
#include <string>

class ModelOptions;
class SimulationInterface;


//! Ramp a variable in quasistationary mode
class TBDLLOCAL Ramp
{

  public:

    //! Constructor
    explicit Ramp(const ModelOptions& options,
        const std::vector<SimulationInterface*>& simulations =
        std::vector<SimulationInterface*>());


    //! Destructor
    ~Ramp(void);

    //! Do the ramp
    void ramp(void);

    //! Ramp to a given goal
    void ramp(double goal);



  private:

    //! The simulations that should be solved
    std::vector<SimulationInterface*> _simulations;

    //! A pointer to the variable
    std::string _variable;

    //! The simulation goal
    double _goal;

    //! The last simulation value
    double _last;

    //! The initial step size 0 < _initial_step < 1
    double _initial_step;

    //! The minimum step size 0 < _min_step < 1
    double _min_step;

    //! The maximum step size 0 < _max_step < 1
    double _max_step;

    //! Whether to plot results or not
    bool _plot_data;


    //! The ids of the remembered solutions
    std::vector<ID> _old_sol_ids;

};


//
// inline methods
//

inline
void
Ramp::ramp(double goal)
{
  _goal = goal;
  ramp();
}




#endif // _SWEEP_H_
