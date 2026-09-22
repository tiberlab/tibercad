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
 * 
 * In case of simplices, we use the standard definition of Whitney
 * forms, which are the same as the FE basis functions for 0-forms,
 * and for 1-forms are defined as
 * \f$\lambda_{ij}=\mathcal{N}_i\mathrm{d\mathcal{N}_j - 
 * \mathcal{N}_j\mathrm{d}\mathcal{N}_i\f$ where 
 * \f$\mathcal{N}_i\f$ is the FE basis function associated to
 * node \f$i\f$.
 * 
 * For non-simplices, this construction in general satisfies interpolation
 * property in the sense that the integral of the 1-form along the
 * edge \f$e_{ij}\f$ is \f$\int_{e_{ij}}\lambda_{kl} = \delta_{ik}\delta_{jl}\f$.
 * However, these 1-forms are inconsistent in the sense that they do not
 * necessarily reproduce constant 1-forms, i.e. the diagram
 * \f$d\circ\mathcal{I}_0 \ne \mathcal{I}_1\circ d_0\f$ does not commute.
 * 
 * For 0-forms, the interpolants are \f$\mathcal{N}_i = \phi_i(x)\f$
 * where \f$\phi_i\f$ is the FEM basis function associated to node
 * \f$i\f$.
 * 
 * In non-simplicial elements, we define 1-form basis functions piecewise
 * by a suitable subdivision of the element into simplices. For example,
 * in a quadrilateral we interpolate in the support formed by the two
 * triangles obtained by diagonal subdivision that share an edge. In
 * this case, the 1-form basis functions are discontinuous across the
 * diagonal of the quad. They are defined as the average of the two
 * Whitney 1-forms on the two triangles, with the diagonal cochain
 * eliminated via Stokes' theorem. This construction is first-order
 * accurate, and is consistent in the sense that it reproduces constant 1-forms. 
 *
 * In this sense, this class provide Whitney-like interpolation forms,
 * which only for simplices are the standard Whitney forms.
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
     * 
     * \param elem the element
     * \param points the points at which to evaluate the interpolants
     * \param local_coordinates true if points are in local coordinates
     */
    void reinit(const libMesh::Elem& elem,
        const std::vector<libMesh::Point>& points, bool local_coordinates = false);

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
     * of the element. Note that 1D elements are assumed to have
     * a single edge, although libMesh doesn't currently assign
     * them one. Also, the 1-forms are returned as RealGradients,
     * containing the coefficients to the coordinate 1-forms dx, dy, dz.
     */
    const std::vector<std::vector<libMesh::RealGradient>>& get_1forms(void) const;

    /*!
     * \brief Retrieve the 1-cells, inclduing virtual ones for non-simplices
     *
     * \return the primal 1-cells, i.e. the edges of the element, as pairs of node indices.
     * This includes the virtual 1-cells due to logical subdivision
     * of non-simplicial elements. It corresponds to the incidence matrix,
     * but is stored as pairs of node indices. The order is the same as in get_1forms().
     */
    const std::vector<std::pair<unsigned int, unsigned int>>& get_1cells(void) const;

    /*!
     * \brief Retrieve the points in real coordinates
     */
    const std::vector<libMesh::Point>& get_xyz(void) const;


  private:

    /*!
     * \brief The current element
     * The pointer is guarantueed to be non-null after
     * reinit() is called.
     */
    const libMesh::Elem* _elem = nullptr;

    /*!
     * \brief the primal 1-cells, i.e. the edges of the element
     * This includes the virtual 1-cells due to logical subdivision
     * of non-simplicial elements. It corresponds to the incidence matrix,
     * but is stored as pairs of node indices.
     * This is needed for the Whitney interpolation of 1-forms, which are defined
     * on the primal 1-cells, and not on the edges of the element, and thus might
     * need evaluation of 1-cochains on additional virtiual 1-cells.
     */
    std::vector<std::pair<unsigned int, unsigned int>> _primal_1cells;

    /*!
     * \brief The 0-forms
     */
    std::vector<std::vector<double>> _w0;

    /*!
     * \brief The 1-forms
     */
    std::vector<std::vector<libMesh::RealGradient>> _w1;

    /*!
     *\brief the points in real coordinates
     */
    std::vector<libMesh::Point> _xyz;


    /*!
     * \brief Setup the primal 1-cells for a given element
     * \param elem the element
     * \param primal_1cells the vector to be filled with the primal 1-cells, as pairs of node indices
     * 
     * \c primal_1cells is filled with the primal 1-cells of the element, which are
     * the edges of the element for simplices, and the edges of the simplicial subelements for non-simplices.
     * The order of the primal 1-cells is consistent with the order of the edges of the element,
     * and is used to define the Whitney interpolation 1-forms. For non-simplices, this means
     * that some primal 1-cells might not correspond to actual edges of the element, but rather
     * to virtual edges that are used to define the Whitney interpolation 1-forms on the simplicial
     * subelements. Real edges precede virtual edges in the in the data structure.
     * This function is called by reinit() to setup the primal 1-cells for the given element.
     */
    void setup_1cells(const libMesh::Elem& elem,
        std::vector<std::pair<unsigned int, unsigned int>>& primal_1cells);

    /*!
     * \brief Find the pair of nodes that form the larger angle
     * \param p0 the first point
     * \param p1 the second point
     * \param p2 the third point
     * \param p3 the fourth point
     * \return the pair of nodes that form the larger angle
     */
    std::pair<unsigned int, unsigned int> larger_angle_pair(const libMesh::Point& p0,
        const libMesh::Point& p1, const libMesh::Point& p2, const libMesh::Point& p3) const;

    /*!
     * \brief Calculate the Whitney forms for a subtriangle
     * \param n1 the first node
     * \param n2 the second node
     * \param n3 the third node
     * \param points the point indices in the subtriangle
     * \param ref_points the reference points
     */
    void calculate_subtriangle_whitney_forms(unsigned int n1, unsigned int n2, unsigned int n3,
        const std::vector<unsigned int>& points, const std::vector<libMesh::Point>& ref_points);


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
const std::vector<libMesh::Point>&
WhitneyInterpolation::get_xyz(void) const
{
  return _xyz;
}

inline
const std::vector<std::pair<unsigned int, unsigned int>>&
WhitneyInterpolation::get_1cells(void) const
{
  return _primal_1cells;
}

#endif // TC_WHITNEYINTERPOLATION_H
