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
 * \file DEC.h
 * \brief Public tiberCAD API header.
 *
 * \note This header is part of the public tiberCAD API.
 *       API version: 3.5
 */



#ifndef TC_DEC_H
#define TC_DEC_H


#include "tibercad/math/WhitneyInterpolation.h"


#include "libmesh/edge.h"



/*!
 * \brief A class to hold Discrete Exterior Calculus related stuff
 * 
 * Thsi class is used to calculate and access quantities needed
 * in DEC, like element center point, incidence matrix, discrete
 * exterior derivative, Hodge star etc.
 *
 * The current implementation uses libMesh's FE order-1 Lagrange
 * element family for the construction of Whiyney interpolation forms.
 *
 * The construction of the dual d-cells can be chosen between 
 * barycentric and circumcentric, or mixed. The default is
 * barycentric, since then the center point is guarantueed to be inside
 * the primal element (d-cell).
 *
 */
class DEC
{

  public:

    /*!
     * \brief Approach for dual element construction
     */
    enum DualConstruction
    {
      BARYCENTRIC,   /*! < barycentric dual construction */
      CIRCUMCENTRIC, /*! < circumcentric dual construction */
      MIXED          /*! < mixed dual construction */
    };
    

    //! Default constructor 
    DEC(void) = delete;

    /*!
     * \brief Constructor with element and scaling
     *
     * \param elem The element to initialize for
     * \param dual_constr The approach for dual element construction
     */
    DEC(const libMesh::Elem& elem,
        DualConstruction dual_constr = BARYCENTRIC);


    /*!
     * \brief initialize the data structures
     */
    void init(void);

    /*!
     * \brief Reinitialize the DEC object for a new element
     *
     * \param elem The element to initialize for
     */
    void reinit(const libMesh::Elem& elem);

    /*! \brief Retrieve the Whitney interpolation object
     * 
     * The internal Whitney forms are calculated at the element center point.
     * 
     * \return A constant reference to the Whitney interpolation object
     */
    const WhitneyInterpolation& get_whitney(void) const { return _whip; }

    //! \brief Retrieve the element center point
    const libMesh::Point& get_center(void) const { return _center; }

    /*!
     * \brief Retrieve the incidence matrix
     * The incidence matrix is the discrete exterior derivative for 0-forms.
     */
    const libMesh::DenseMatrix<double>& get_incidence_matrix(void) const { return _incidence; }

    //! \brief Retrieve the dual volumes
    const std::vector<double>& get_dual_volumes(void) const { return _dual_volumes; }

    //! \brief Retrieve the primal 1-cell vectors
    const std::vector<libMesh::RealGradient>& get_primal_1cells(void) const { return _primal; }

    //! \brief Retrieve the midpoints of the edges
    const std::vector<libMesh::Point>& get_midpoints(void) const { return _midpoints; }

    /*!
     * \brief Get the Hodge star, possibly including a metric factor
      * \param hodge The Hodge star matrix to be filled
      * \param metric The metric matrix to be used, if any
      *
     */
    void get_hodge(libMesh::DenseMatrix<double>& hodge,
        const libMesh::RealTensor& metric = libMesh::RealTensor(1, 0, 0, 0, 1, 0, 0, 0, 1)) const;


    //! Set the dual construction approach
    void set_dual_construction(DualConstruction dual_constr) { _dual_constr = dual_constr; }


    //! Get the incidence as pairs of nodes, in the same order as the incidence matrix
    void get_incidence_pairs(std::vector<std::pair<unsigned int, unsigned int>>& inc) const;


  private:

    //! The Whitney interpolation object
    WhitneyInterpolation _whip;

    /*!
     * \brief The element we are currently working on
     * The pointer is guarantueed to be non-null after
     * the constructor or reinit() is called.
     */
    const libMesh::Elem* _elem = nullptr;

    //! The dual construction approach
    DualConstruction _dual_constr = BARYCENTRIC;

    //! The element center point
    libMesh::Point _center;

    //! The incidence matrix = d0, the discrete exterior derivative for 0-forms
    libMesh::DenseMatrix<double> _incidence;

    //! The dual 0-cell volume contributions for each primal node
    std::vector<double> _dual_volumes;

    //! The primal 1-cell vectors for each edge, ordered as the edges of the element
    std::vector<libMesh::RealGradient> _primal;

    //! The midpoints of the edges, ordered as the edges of the element
    std::vector<libMesh::Point> _midpoints;

    //! \brief Calculate the circumcenter of a given element
    /*!
     * \param elem The element to calculate the circumcenter for
     * \param s The side index to use for the circumcenter calculation, if applicable
     * \return The circumcenter point of the element
     *
     * This function calculates the circumcenter of the given element. If a side index
     * is provided, the circumcenter is calculated such that the points of that side
     * lie on the circumcircle.
     * Circumcenter is not implemented for all element types, and the function may fall
     * back to the barycenter
     */
    libMesh::Point circumcenter(const libMesh::Elem& elem, int s = -1) const;

    void compute_quad_hodge(const libMesh::Elem& elem,
                            libMesh::DenseMatrix<double>& H,
                            const libMesh::RealTensor& metric) const;

};


#endif // TC_DEC_H
