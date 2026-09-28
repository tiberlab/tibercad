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



#include "libmesh/point.h"
#include "libmesh/vector_value.h"

#include <vector>
#include <utility>




//! Whitney interpolation forms for DEC
/*!
 * Here we define Whitney interpolation forms, or more precisely
 * the respective coefficient functions, that are needed in
 * the Discrete Exterior Calculus (DEC) formulation of PDEs.
 * Primal interpolation of 0- and 1-forms is implemented for simplices,
 * i.e. an edge (2 nodes), a triangle (3 nodes) or a tetrahedron (4 nodes),
 * given directly by their node coordinates.
 *
 * The standard definition of Whitney forms is used, which for 0-forms
 * are the barycentric coordinate functions \f$\mathcal{N}_i\f$ of the
 * simplex, and for 1-forms are defined as
 * \f$\lambda_{ij}=\mathcal{N}_i\mathrm{d}\mathcal{N}_j -
 * \mathcal{N}_j\mathrm{d}\mathcal{N}_i\f$, where
 * \f$\mathcal{N}_i\f$ is the barycentric coordinate function associated
 * to node \f$i\f$.
 *
 * This class does not depend on libMesh's \c Elem class or mesh
 * connectivity: it only needs the node coordinates of the simplex to
 * interpolate on. Non-simplicial elements (e.g. quadrilaterals or
 * hexahedra) are not handled here; logically subdividing them into
 * sub-simplices, evaluating this class on each sub-simplex, and
 * combining the results is the responsibility of the DEC class.
 */
class WhitneyInterpolation
{

  public:

    /*!
     * \brief Default constructor
     */
    WhitneyInterpolation(void) = default;

    /*!
     * \brief Recalculate interpolants for a given simplex and points
     *
     * \param nodes the nodes of the simplex, in real coordinates: 2 nodes
     *        for an edge (1-simplex), 3 for a triangle (2-simplex), or
     *        4 for a tetrahedron (3-simplex)
     * \param points the points at which to evaluate the interpolants, in
     *        real coordinates
     */
    void reinit(const std::vector<libMesh::Point>& nodes,
        const std::vector<libMesh::Point>& points);

    /*!
     * \brief Retrieve the 0-forms
     *
     * The first vector index refers to the interpolant, ordered as the
     * nodes passed to reinit(); the second to the point.
     */
    const std::vector<std::vector<double>>& get_0forms(void) const;

    /*!
     * \brief Retrieve the 1-forms
     *
     * The first vector index refers to the interpolant, ordered as
     * returned by get_1cells(); the second to the point. The 1-forms are
     * returned as RealGradients, containing the coefficients to the
     * coordinate 1-forms dx, dy, dz.
     */
    const std::vector<std::vector<libMesh::RealGradient>>& get_1forms(void) const;

    /*!
     * \brief Retrieve the 1-cells of the simplex
     *
     * \return the 1-cells, i.e. the edges of the simplex, as pairs of
     * node indices (into the \c nodes array passed to reinit()). The
     * edges are listed as all pairs (i,j) with i<j, in lexicographic
     * order. The order is the same as in get_1forms().
     */
    const std::vector<std::pair<unsigned int, unsigned int>>& get_1cells(void) const;


  private:

    /*!
     * \brief The simplex nodes, in real coordinates
     */
    std::vector<libMesh::Point> _nodes;

    /*!
     * \brief The gradients of the barycentric coordinate functions
     * (0-forms), one per node. These are constant over the simplex.
     */
    std::vector<libMesh::RealGradient> _grad;

    /*!
     * \brief the 1-cells, i.e. the edges of the simplex, as pairs of
     * node indices, in lexicographic order
     */
    std::vector<std::pair<unsigned int, unsigned int>> _cells;

    /*!
     * \brief The 0-forms
     */
    std::vector<std::vector<double>> _w0;

    /*!
     * \brief The 1-forms
     */
    std::vector<std::vector<libMesh::RealGradient>> _w1;


    /*!
     * \brief Setup the 1-cells for a simplex with the given number of nodes
     */
    void setup_1cells(void);

    /*!
     * \brief Compute the (constant) gradients of the barycentric
     * coordinate functions of the simplex
     */
    void compute_gradients(void);

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


inline
const std::vector<std::pair<unsigned int, unsigned int>>&
WhitneyInterpolation::get_1cells(void) const
{
  return _cells;
}

#endif // TC_WHITNEYINTERPOLATION_H
