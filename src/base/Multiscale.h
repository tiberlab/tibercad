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
 * \file Multiscale.h
 * \brief Private tiberCAD header.
 *
 * \internal
 */


#ifndef TC_MULTISCALE_H
#define TC_MULTISCALE_H


#include "tibercad/base/ModelOptions.h"
#include "tibercad/base/TypeDefs.h"

#include <set>

class SimulationInterface;

//! Class for the handling of multiscale couplings
class Multiscale
{

  public:

    //! The coupling method
    enum Method
    {
      UNKNOWN, //!< invalid method
      OVERLAP, //!< overlap method
      BRIDGE   //!< bridge method
    };

    //! Constructor
    explicit Multiscale(const ModelOptions& options);

    //! Reinit the multiscale coupling
    /*!
     * Can be called before solving \c sim.
     * If \c sim is one of the multiscale coupled models,
     * actions will be taken according to multiscale method.
     */
    void reinit(SimulationInterface* sim);


  private:

    //! The options
    ModelOptions _options;

    //! The coupling method
    Method _method;

    //! The macroscale model
    SimulationInterface* _macromodel;

    //! The microscale model
    SimulationInterface* _micromodel;

    //! The restriction of the macro domain
    std::set<ID> _restricted_macro_domain;

    //! The restricted variables
    std::vector<std::string> _restricted_variables;

    //! Initialize
    void init(void);

};


#endif // TC_MULTISCALE_H
