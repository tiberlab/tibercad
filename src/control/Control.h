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
 * \file Control.h
 * \brief Private tiberCAD header.
 *
 * \internal
 */


#ifndef _CONTROL_H_
#define _CONTROL_H_

#include "tibercad/base/TypeDefs.h"
#include "tibercad/base/IDSet.h"
#include "tibercad/base/tiber_dll.h"

#include <map>
#include <list>
#include <vector>
#include <set>
#include <iostream>
#include <cassert>

class Device;
class Database;
class ModelOptions;
class SimulationInterface;

//! The control module of TiberCAD
/*!
 * This class is responsible for the setup of the device and
 * the flow control of the simulation
 */
class TBDLLOCAL Control
{

  public:


    //! The constructor
    Control(void);


    //! The destructor
    ~Control(void);


    //! Set the input file
    void set_inputfile(const std::string& inputfile);


    //! Initialize the device and the simulations
    /*!
     * This method calls create_device(), create_materials() and
     * setup_models()
     */
    void init(void);


    //! Runs the simulation
    /*!
     * Calls solve() of the simulations according to the rules given
     * in the input file
     */
    void run_simulation(void);


    //! Get the simulation time coordinate
    double get_time(void);


    //! Plots the results of all simulation models
    void plot_all(void);


    //! Get a reference to the device
    Device& get_device(void);



    //! Get the directory where to put output files
    const std::string& get_output_dir(void) const;



    //! Get the output format
    /*!
     * \return a string that identifies the type ouf output files
     * to generate
     *
     * Currently the following formats are supported:
     * \li \c gmv for GMV
     * \li \c ise for Tecplot
     * \li \c gnu for GnuPlot
     */
    const std::string& get_output_format(void) const;



  private:

    /*
    //! A signal handler
    class SignalHandler
    {

      public:

        typedef struct sigaction SigAction;

        static void set_control(Control* ctrl);

        static void activate_sigint(void);
        static void deactivate_sigint(void);

      private:

        SignalHandler(void);

        static Control* _ctrl;

        static SigAction _old_int_action;

        static void sigint(int sig);
    };
    */


    //! The inputfile
    std::string _inputfile;


    //! The device we control
    Device* _device = nullptr;


    //! The time variable
    double _time = 0.0;


    //! The database we use
    Database* _database = nullptr;


    //! The list of simulations to be solved
    std::vector<std::string> _solve_list;


    //! The directory where to put output
    std::string _outputdir = ".";


    //! The output format
    /*!
     * see get_output_format() for a detailed description
     */
    std::string _output_format;


    //! Setup global simulation options
    void setup_globals(const ModelOptions& opts);


    //! Setup a model
    void setup_module(Device* device, const ModelOptions& opts);


};


//
// inline members
//


inline
void
Control::set_inputfile(const std::string& inputfile)
{
  _inputfile = inputfile;
}


inline
Device&
Control::get_device(void)
{
  return *_device;
}



inline
const std::string&
Control::get_output_dir(void) const
{
  return _outputdir;
}




inline
const std::string&
Control::get_output_format(void) const
{
  return _output_format;
}





#endif // _CONTROL_H_
