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
 * \file Sweep.h
 * \brief Private tiberCAD header.
 *
 * \internal
 */


#ifndef TC_SWEEP_H
#define TC_SWEEP_H

#include "tibercad/module/SimulationInterface.h"

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
class TC_DLLOCAL Sweep : public SimulationInterface
{

  public:
    
    //! Destructor
    virtual ~Sweep(void);

    //! Create a Sweep object
    static Sweep* create(const ModelOptions& options);


    //! Get the list of inner simulations
    const std::vector<SimulationInterface*>&
      get_simulations(void) const;



  protected:

    //! Constructor
    Sweep(const ModelOptions& options);

    /*! \copydoc SimulationInterface::do_init() */
    virtual void do_init(void);


    // / *! \copydoc SimulationInterface::do_equilibrium() */
    //virtual void do_equilibrium(void);

    
    /*! \copydoc SimulationInterface::do_solve()
     *
     * solve() of the associated SimulationInterface objects will
     * be called repeatedly for all sweep values
     */
    virtual void do_solve(void);
   
    //! Plot final results.
    virtual void do_plot(void);


    //! Abused here for the inheritance of variables
    virtual void do_print_info(void);


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

    //! Types of sweeps
    enum Type
    {
      LINEAR,
      LOG
    };

    //! The simulations for which we do the sweep
    std::vector<SimulationInterface*> _simulations;

    //! A pointer to the sweepable variable
    /*!
     * We use a vector of string so that we can sweep
     * multiple variables together, if needed.
     */
    std::vector<std::string> _variable;

    //! The type of the sweep
    Type _type {LINEAR};

    //! A vector containing all the simulation values
    std::vector<double> _values;

    //! The current simulation goal
    double _goal {0.0};


    //! When to write results to file
    std::vector<bool> _plot_data;

    //! Whether to append sweep name to filename suffix
    bool _append_sweep_name_to_datafile_name {false};

    //! The current filename suffix
    std::string _suffix {""};


    //! The ids of the remembered solutions
    std::map<ID, std::vector<ID> > _remembered_sol_ids;

    //! The remembered sweep value
    double _remembered_sweep_value {0.0};


    //! Ignore solver failures
    bool _ignore_failures {false};


    //! Remember the current solution
    ID remember_solution(void);


    //! Do the sweep
    void do_sweep(std::vector<double>& values,
        std::vector<std::ofstream*>& plotfiles,
        std::vector<std::map<double, std::vector<double> > >& sweep_data);


    //! Write the global data to file
    void write_global_data(SimulationInterface& simulation, std::ofstream*& file);


};




#endif // TC_SWEEP_H
