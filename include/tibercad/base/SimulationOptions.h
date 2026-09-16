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
 * \file SimulationOptions.h
 * \brief Public tiberCAD API header.
 *
 * \note This header is part of the public tiberCAD API.
 *       API version: 3.5
 */


#ifndef TC_SIMULATIONOPTIONS_H
#define TC_SIMULATIONOPTIONS_H

#include "tibercad/base/tiber_dll.h"
#include <string>

class ModelOptions;

//! Common options for set of simulations
class SimulationOptions
{

  public:


    //! The ambient temperature in K
    static double temperature;
    static double& temp;
    static double& T;

    //! Whether or not to consider incomplete ionization
    static bool incomplete_ionization;


    //! Tell verbosity
    static int verbose(void);


    //! Initialize the simulation options
    static void initialize(const ModelOptions& opts) TC_DLLOCAL;

    //! A path for temporary data
    static std::string scratch_path;

  private:

    SimulationOptions(void) TC_DLLOCAL {};
    
    //! Level of verbosity
    static int _verbose;

};



#endif //_SIMULATIONOPTIONS_H_
