/*  
 * This file is part of the tiberCAD module pvmodule.
 *
 * tiberCAD modules are licensed under the GNU General Public License v3.
 *
 * tiberCAD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published
 * by the Free Software Foundation, either version 3 of the License,
 * or (at your option) any later version.
 *
 * tiberCAD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with tiberCAD. If not, see <http://www.gnu.org/licenses/>.
 */

/*!
 * \file DegradationActiveFraction.h
 * \brief tiberCAD pvmodule module header.
 *
 * \note This file is part of module pvmodule.
 */

#ifndef TC_DEGRADATIONACTIVEFRACTION_H
#define TC_DEGRADATIONACTIVEFRACTION_H

#include "DegradationModel.h"
#include "tibercad/module/SolutionProvider.h"

/*!
 * \brief An example for degradation and its impact on performance
 *
 * This class implements an equiv. circuit parameter dependency 
 * on the fraction of active (non-degraded) perovskite, calculated from a kinetic model.
 * 
 * J = J0 x f^n
 * Rs = Rs0 / f^n
 * 
 * Here f is the fraction of active (non-degraded) perovskite, which has to be provided by
 * another module.
 */
class DegradationActiveFraction : public DegradationModel
{

  public:

    virtual ~DegradationActiveFraction(void) = default;


  protected:

    //! Private constructor
    DegradationActiveFraction(const ModelOptions& options);

    virtual void do_init(void) override;

    virtual void do_degrade_params(const libMesh::Elem* elem,
                                   const libMesh::Point& p,
                                   DegradationModel::Parameters& params) const final;


  private:

    //! The exponent in the photocurrent degradation dependency
    double _exponent_ph = 1;

    //! The exponent in the series resistance degradation dependency
    double _exponent_rs = 1;

    //! From where to get active perovskite fraction
    SolutionProvider _kinetic_model;

};


#endif // TC_DEGRADATIONACTIVEFRACTION_H
