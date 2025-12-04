/*
 * This file is part of tiberCAD.
 *
 * tiberCAD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * tiberCAD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with tiberCAD. If not, see <https://www.gnu.org/licenses/>.
 */

/*!
 * \file Ramp.h
 * \brief Private tiberCAD header.
 *
 * \internal
 */


#ifndef _RAMP_H_
#define _RAMP_H_

#include "tibercad/base/tiber_dll.h"
#include "tibercad/base/TypeDefs.h"

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

    //! Get the current value



  private:

    //! The simulations that should be solved
    std::vector<SimulationInterface*> _simulations;

    //! A pointer to the variable
    /*!
     * We use a vector of string so that we can sweep
     * multiple variables together, if needed.
     */
    std::vector<std::string> _variable {""};

    //! The simulation goal
    double _goal {0.0};

    //! The last simulation value
    double _last {0.0};

    //! The initial relative step size 0 < _initial_step < 1
    double _initial_step {1};

    //! The minimum relative step size 0 < _min_step < 1
    double _min_step {0.01};

    //! The maximum relative step size 0 < _max_step < 1
    double _max_step {1};


    //! The initial absolute step size 0 < _initial_step < 1
    double _initial_abs_step {std::numeric_limits<double>::max()};

    //! The minimum absolute step size 0 < _min_step < 1
    double _min_abs_step {std::numeric_limits<double>::min()};

    //! The maximum absolute step size 0 < _max_step < 1
    double _max_abs_step {std::numeric_limits<double>::max()};


    //! Whether to plot results or not
    bool _plot_data {false};


    //! The ids of the remembered solutions
    std::vector<ID> _old_sol_ids;

};





#endif // _SWEEP_H_
