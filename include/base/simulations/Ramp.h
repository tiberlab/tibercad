// $Id$

#ifndef _RAMP_H_
#define _RAMP_H_

#include "TypeDefs.h"

#include <set>
#include <vector>
#include <string>

class ModelOptions;
class SimulationInterface;


//! Ramp a variable in quasistationary mode
class Ramp
{

  public:

    //! Constructor
    explicit Ramp(const ModelOptions& options);

    //! Destructor
    ~Ramp(void);

    //! Do the ramp
    void ramp(void);


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


    //! The ids of the remembered solutions
    std::map<ID, std::vector<ID> > _remembered_sol_ids;

};


//
// inline methods
//




#endif // _SWEEP_H_
