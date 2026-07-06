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
 * \file WhitneyInterpolation.h
 * \brief Public tiberCAD API header.
 *
 * \note This header is part of the public tiberCAD API.
 *       API version: 3.5
 */



#ifndef TC_WHITNEYINTERPOLATION_H
#define TC_WHITNEYINTERPOLATION_H



#include "libmesh/fe_base.h"

#include <vector>
#include <memory>




//! Whitney interpolation forms for DEC
/*!
 * Here we define Whitney interpolation forms, or more precisely
 * the respective coefficient functions, that are needed in
 * the Discrete Exterior Calculus (DEC) formulation of PDEs.
 * For now, primal interpolation of 0- and 1-forms is implemented,
 * which is based on the Finite Element Basis functions.
 * For 0-forms, the interpolants are \f$\mathcal{N}_i = \phi_i(x)\f$
 * where \f$\phi_i\f$ is the FEM basis function associated to node
 * \f$i\f$.
 * For 1-forms, they are \f$\lambda_{ij}=\mathcal{N}_i\mathrm{d}\mathcal{N}_j
 * - \mathcal{N}_j\mathrm{d}\mathcal{N}_i\f$.
 *
 */
class WhitneyInterpolation
{

  public:

    /*!
     * \brief Default constructor
     */
    WhitneyInterpolation(void) = default;
    
    /*!
     * \brief Recalculate interpolants for given element and points
     */
    void reinit(const libMesh::Elem& elem,
        const std::vector<libMesh::Point>& points);

    /*!
     * \brief Retrieve the 0-forms
     *
     * The first vector index refers to the point, the second
     * to the interpolant. The latter are ordered as the nodes
     * of the element.
     */
    const std::vector<std::vector<double>>& get_0forms(void) const;

    /*!
     * \brief Retrieve the 1-forms
     *
     * The first vector index refers to the point, the second
     * to the interpolant. The latter are ordered as the edges
     * of the element. Not that 1D elements are assumed to have
     * a single edge, although libMesh doesn't currently assign
     * them one. Also, the 1-forms are returned as RealGradients,
     * containing the coefficients to the coordinate 1-forms dx, dy, dz.
     */
    const std::vector<std::vector<libMesh::RealGradient>>& get_1forms(void) const;



  private:

    /*!
     * \brief The 0-forms
     */
    std::vector<std::vector<double>> _w0;

    /*!
     * \brief The 1-forms
     */
    std::vector<std::vector<libMesh::RealGradient>> _w1;
};


inline
const std::vector<std::vector<double>>&
WhitneyInterpolation::get_0forms(void) const
{
  return _w0; 
}


inline
const std::vector<std::vector<libMesh::RealGradient>>&
WhitneyInterpolation::get_1forms(void) const
{
  return _w1;
}

#endif // TC_WHITNEYINTERPOLATION_H
