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
 * \file SolutionProvider.h
 * \brief Public tiberCAD API header.
 *
 * \note This header is part of the public tiberCAD API.
 *       API version: 3.5
 */


#ifndef _SOLUTIONPROVIDER_H_
#define _SOLUTIONPROVIDER_H_

#include "tibercad/base/TypeDefs.h"

#include <utility>

class SimulationInterface;

/*!
 * \brief A class providing access to solutions in other modules
 */
class SolutionProvider : public std::pair<SimulationInterface*, ID>
{
  public:

    SolutionProvider()
      : std::pair<SimulationInterface*, ID>(nullptr, INVALID_ID)
    {}


    SolutionProvider(SimulationInterface* si, ID id)
      : std::pair<SimulationInterface*, ID>(si, id)
    {}


    explicit SolutionProvider(const std::pair<SimulationInterface*, ID>& in)
      : std::pair<SimulationInterface*, ID>(in)
    {}


    //! Get the pointer to the simulation object
    SimulationInterface* simulation(void) const
    {
      return(this->first);
    }

    //! Get the ID
    ID id(void) const
    {
      return(this->second);
    }

    //! Check if it is valid
    bool is_valid(void) const
    {
      return((this->first != nullptr) && (this->second != INVALID_ID));
    }

};


#endif // _SOLUTIONPROVIDER_H_
