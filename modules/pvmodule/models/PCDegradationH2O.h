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
 * \file PCDegradationH2O.h
 * \brief tiberCAD pvmodule module header.
 *
 * \note This file is part of module pvmodule.
 */

#ifndef TC_PCDEGRADATIONH2O_H
#define TC_PCDEGRADATIONH2O_H

#include "Photocurrent.h"
#include "tibercad/module/SolutionProvider.h"

/*!
 * \brief An example for photocurrent dependency on H20 concentration
 *
 * This class implements a photocurrent dependency based on the
 * assumption, that presence of water in a perovskite cell leads
 * to material degradation and reduction of optically active
 * material.
 */
class PCDegradationH2O : public Photocurrent
{

  public:

    ~PCDegradationH2O(void) {};

    static PCDegradationH2O* create(const ModelOptions& options);


  protected:

    PCDegradationH2O(const ModelOptions& options);

    virtual void do_init(void) final;

    virtual double do_get_photocurrent(const libMesh::Elem* elem,
                                       const libMesh::Point& p) const final;


  private:

    //! The current degradation factor
    /*!
     * The formula is based on a fit of the data in [ref].
     * It uses a generalized logistic function.
     */
    double degradation_factor(double humidity) const;

    //! The initial (undegraded) photocurrent
    double _initial_current = 0.02;

    //! The reference humidity in the degradation fit
    double _RH_ref = 72.1;

    //! The exponent in the degradation fit
    double _exponent = 8.28;

    //! From where to get relative humidity
    SolutionProvider _humidity_model;

};


#endif // TC_PCDEGRADATIONH2O_H
