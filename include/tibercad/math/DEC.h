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
#include "libmesh/dense_matrix.h"
#include "libmesh/tensor_value.h"

#include <map>



/*!
 * \brief A class to hold Discrete Exterior Calculus related stuff
 * 
 * This class is used to calculate and access quantities needed
 * in DEC, like element center point, incidence matrix, discrete
 * exterior derivative, Hodge star etc.
 *
 * For simplices, WhitneyInterpolation is used directly on the element's
 * own nodes to construct the standard Whitney interpolation forms.
 * The construction of the dual d-cells can be chosen between
 * barycentric and circumcentric, or mixed. The default is
 * barycentric, since then the center point is guarantueed to be inside
 * the primal element (d-cell).
 * 
 * For non-simplices, this class logically subdivides the element into
 * sub-simplices (a quadrilateral into two triangles, a pyramid into two
 * tetrahedra, a prism into three tetrahedra, a hexahedron into six
 * tetrahedra; see setup_subdivision()), and uses WhitneyInterpolation on
 * the sub-simplex containing a given point to construct piecewise Whitney
 * interpolation
 * 1-forms on the non-simplicial element; these are discontinuous across
 * the internal (virtual) sub-simplex boundaries. The 0-forms, on the
 * other hand, use the element's own native (e.g. bilinear/trilinear)
 * Lagrange basis functions, which are smooth over the whole element.
 *
 * Alternatively, the Hodge star can be calculated using a mimetic finite
 * difference (MFD) approach, which is consistent with linear fields, graph
 * compatible and might have better numerical properties, at a cost of
 * inconsistency with the interpolation forms.
 *
 */
class DEC
{

  public:

    /*!
     * \brief Approach for dual element construction on simplices
     */
    enum DualConstruction : unsigned int
    {
      BARYCENTRIC,   /*! < barycentric dual construction */
      CIRCUMCENTRIC, /*! < circumcentric dual construction */
      MIXED          /*! < mixed dual construction */
    };

    /*!
     * \brief Approach for Hodge star construction
     */
    enum HodgeConstruction : unsigned int
    {
      INTERPOLATION, /*! < Hodge star construction using Whitney interpolation */
      MFD            /*! < Hodge star construction using mimetic finite differences */
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
        DualConstruction dual_constr = BARYCENTRIC,
        HodgeConstruction hodge_constr = INTERPOLATION);


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

    /*!
     * \brief Evaluate the 0- and 1-form interpolants at the given points
     *
     * The 0-forms are the element's own native (e.g. bilinear/trilinear)
     * Lagrange basis functions, smooth over the whole element. For the
     * 1-forms, each point is assigned to the sub-simplex of the (possibly
     * logically subdivided, see get_primal_1cells()) element that
     * contains it, and the piecewise Whitney interpolants are evaluated
     * there. The results are accessible via get_0forms() and
     * get_1forms(), indexed in the global node / primal 1-cell numbering
     * of the element.
     *
     * \param points the points at which to evaluate the interpolants
     * \param reference_coord if true, \c points are given in reference
     *        (local) coordinates on the element, and are mapped to real
     *        coordinates via libMesh::FEMap before evaluation; otherwise
     *        \c points are already in real coordinates. Either way, the
     *        real coordinates used are available via get_xyz().
     */
    void reinit_forms(const std::vector<libMesh::Point>& points,
                      bool reference_coord = false);

    /*!
     * \brief Retrieve the points (in real coordinates) used by the last
     * call to reinit_forms()
     */
    const std::vector<libMesh::Point>& get_xyz(void) const { return _xyz; }

    /*!
     * \brief Retrieve the 0-forms computed by the last call to reinit_forms()
     *
     * The first vector index refers to the interpolant, ordered as the
     * nodes of the element; the second to the point.
     */
    const std::vector<std::vector<double>>& get_0forms(void) const { return _w0; }

    /*!
     * \brief Retrieve the 1-forms computed by the last call to reinit_forms()
     *
     * The first vector index refers to the interpolant, ordered as
     * get_primal_1cells() (including virtual 1-cells for logically
     * subdivided, non-simplicial elements); the second to the point.
     */
    const std::vector<std::vector<libMesh::RealGradient>>& get_1forms(void) const { return _w1; }

    /*!
     * \brief Retrieve the element center point
     *
     * On simplices, the center point corresponds to the barycenter or circumcenter,
     * depending on the dual construction approach. On non-simplices, the center
     * point might depend on the simplicial subdivision.
     * 
     * \return A constant reference to the element center point
     */
    const libMesh::Point& get_center(void) const { return _center; }

    /*!
     * \brief Retrieve the incidence matrix
     *
     * The incidence matrix is the discrete exterior derivative for 0-forms.
     * Note that the incidence matrix might contain entries for internal edges
     * of the element, which are not part of its primal 1-cells. This is the case
     * for non-simplicial elements, where the primal 1-cells are constructed from
     * the edges of the simplicial subelements.
     * 
     * \return A constant reference to the incidence matrix
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
        const libMesh::RealTensor& metric = libMesh::RealTensor(1, 0, 0, 0, 1, 0, 0, 0, 1));


    //! Set the dual construction approach
    void set_dual_construction(DualConstruction dual_constr) { _dual_constr = dual_constr; }

    //! Set the Hodge construction approach
    void set_hodge_construction(HodgeConstruction hodge_constr) { _hodge_constr = hodge_constr; } 


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

    //! The Hodge construction approach
    HodgeConstruction _hodge_constr = INTERPOLATION;

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

    /*!
     * \brief The primal 1-cells, i.e. the edges of the element, as pairs
     * of (global) node indices
     *
     * For a simplex, this simply lists all edges of the element. For a
     * logically subdivided (non-simplicial) element, real edges of the
     * element precede the virtual edges introduced by the subdivision.
     */
    std::vector<std::pair<unsigned int, unsigned int>> _primal_1cells;

    /*!
     * \brief The sub-simplices the element is logically subdivided into
     *
     * Each entry lists the (global) node indices of one sub-simplex. For
     * a simplex element, there is a single entry listing all its nodes.
     */
    std::vector<std::vector<unsigned int>> _sub_simplices;

    //! Lookup from a (sorted) node pair to its row in _primal_1cells
    std::map<std::pair<unsigned int, unsigned int>, unsigned int> _cell_index;

    //! The 0-forms computed by the last call to reinit_forms()
    std::vector<std::vector<double>> _w0;

    //! The 1-forms computed by the last call to reinit_forms()
    std::vector<std::vector<libMesh::RealGradient>> _w1;

    //! The points (in real coordinates) used by the last call to reinit_forms()
    std::vector<libMesh::Point> _xyz;

    /*!
     * \brief Calculate the center point of a given element 
     *
     * \param elem The element to calculate the center for
     * \return The center point of the element
     * This function calculates the center point of the given element. 
     * The center point is calculated according to the chosen approach.
     * The default is the barycenter.
     */
    libMesh::Point get_center(const libMesh::Elem& elem) const;

    /*!
     * \brief Calculate the circumcenter of a given element
     *
     * \param elem The element to calculate the circumcenter for
     * \return The circumcenter point of the element
     *
     * This function calculates the circumcenter of the given element. 
     * Circumcenter is not implemented for all element types, and the function may fall
     * back to the barycenter
     */
    libMesh::Point circumcenter(const libMesh::Elem& elem) const;

    /*!
     * \brief Compute intersection of quad diagonals
     * \param elem The quadrilateral element
     * \return The intersection point of the diagonals
     */
    libMesh::Point diagonal_intersection(const libMesh::Elem& elem) const;

    /*!
     * \brief Setup the logical subdivision of an element into sub-simplices
     *
     * For a simplex, \c _sub_simplices contains a single entry listing all
     * of its nodes, and no virtual 1-cells are introduced.
     *
     * For a quadrilateral, the element is split into two triangles by the
     * diagonal connecting the pair of nodes with the larger sum of
     * subtended angles (see larger_angle_pair()).
     *
     * For a pyramid (a quadrilateral base, nodes 0-3, plus an apex, node
     * 4), the element is split into two tetrahedra by the same diagonal
     * of the base as for a quadrilateral.
     *
     * For a triangular prism (a bottom face 0-1-2 and a top face 3-4-5,
     * with real vertical edges i-(i+3)), the element is split into three
     * tetrahedra, pivoting on node 0.
     *
     * For a hexahedron, the element is split into six tetrahedra sharing
     * the main diagonal between nodes 0 and 6.
     */
    void setup_subdivision(void);

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
     * \brief Compute the Hodge star using a mimetic finite difference approach
     * \param elem The element
     * \param H The Hodge star matrix to be filled
     * \param metric The metric tensor to be used in the computation
     *
     * This function computes the Hodge star for 2D and 3D elements using a
     * mimetic finite difference approach.
     * It takes into account the provided metric tensor and fills the Hodge star
     * matrix accordingly.
    */
    void compute_hodge_mfd(const libMesh::Elem& elem,
                           libMesh::DenseMatrix<double>& H,
                           const libMesh::RealTensor& metric) const;
};


#endif // TC_DEC_H
