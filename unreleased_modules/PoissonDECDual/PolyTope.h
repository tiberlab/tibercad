#ifndef _POLYTOPE_H_
#define _POLYTOPE_H_

#include "libmesh/point.h"
#include "libmesh/vector_value.h"

#include <vector>


/*!
 * A polytope defined by its nodes. The main purpose
 * of the class is to provide interpolation methods
 * on general polytopes for use in 2D discrete exterior
 * calculus discretizations. It is assumed that the nodes
 * are ordered counter-clock wise.
 * 
 * \brief A polytope
 * \author Matthias Auf der Maur
 * \date 2025-03-20
 */
class PolyTope
{
  public:

    PolyTope(void) = default;
    PolyTope(PolyTope&&) = default;
    PolyTope(const PolyTope&) = default;
    PolyTope& operator=(PolyTope&&) = default;
    PolyTope& operator=(const PolyTope&) = default;
    ~PolyTope(void) = default;

    //! Clear all nodes
    void clear(void);

    //! Reserve memory for \c n nodes
    void reserve(unsigned int n);

    /*!
     * Add nodes to the polytope object. Nodes
     * have to be added in correct order such as
     * to form a counter-clock wise oriented polytope.
     *
     * \param point libMesh::Point to be added
     * 
     * \brief Add a node
     */
    void add_node(const libMesh::Point& point);

    //! Get access to node list
    const std::vector<libMesh::Point>& get_nodes(void) const;

    /*!
     * Get the barycentric interpolation and its
     * derivatives at \c x . From these, we can construct the
     * necessary 0- and 1-Whitney forms.
     *
     * \param x the coordinates, assumed inside the Polytope
     * \param grad the gradient
     * 
     * \brief Get the interpolation and its derivatives
     */
    double get_interpolation(const libMesh::Point& x,
                             libMesh::RealGradient& grad) const;

  private:

    //! The ordered list of nodes
    std::vector<libMesh::Point> _nodes;

};

#endif // _POLYTOPE_H_
