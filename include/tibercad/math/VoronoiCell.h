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
 * \file VoronoiCell.h
 * \brief Public tiberCAD API header.
 *
 * \note This header is part of the public tiberCAD API.
 *       API version: 3.5
 */


#ifndef TC_VORONOICELL_H
#define TC_VORONOICELL_H

#include <vector>

// forward declaration of elem
namespace libMesh { class Elem; }


/*!
 * \brief Calculate Voronoi cell quantities on an element
 *
 * Calculates the geometric quantities of the partial Voronoi cells
 * associated to the given element, for all of its nodes. 
 */
class VoronoiCell
{

  public:

    //! A typedef for convenience
    typedef std::pair<unsigned int, unsigned int> NodePair;

    //! Constructor
    explicit VoronoiCell(const libMesh::Elem* elem, double scaling = 1.0);

    //! Constructor
    explicit VoronoiCell(const libMesh::Elem& elem, double scaling = 1.0);

    //! Destructor
    ~VoronoiCell(void);

    //! Calculate the cell quantities
    void calculate(void);

    //! Get the number of edges
    unsigned int n_edges(void) const;

    //! Get the edge indices for a given element node
    //const std::vector<unsigned int>& get_edges(unsigned int n) const;

    //! Get the edge nodes
    const NodePair& get_nodes(unsigned int e) const;

    //! Get the edge length
    double get_edge_length(unsigned int e) const;

    //! Get the edge lengths
    const std::vector<double>& get_edge_lengths(void) const;

    //! Get the area of the cell face associated with given edge
    double get_cell_face_area(unsigned int e) const;

    //! Get the areas of the cell faces
    const std::vector<double>& get_cell_face_areas(void) const;

    //! Get the volume associated to a given node
    double get_volume(unsigned int n) const;

    //! Get the partial cell volumes
    const std::vector<double>& get_volumes(void) const;


  private:

    //! The element this object is related to
    const libMesh::Elem* _elem;

    //! The length scaling
    double _scaling;

    //! The edges, normal to cell faces
    /*!
     * The two numbers in the std::pair are the two nodes,
     * in local indexing.
     */
    std::vector<NodePair> _edges;

    //! The areas of the edges
    std::vector<double> _edge_areas;

    //! The edge lengths
    std::vector<double> _edge_lengths;

    //! The partial cell volumes
    std::vector<double> _volumes;

    //! The association between nodes and edges
    //std::vector<std::vector<unsigned int>> _nodes_to_edges;

    //! Clear all data structures
    void clear(void);

    //! Solve a 2x2 system
    std::pair<double, double> solve_2x2(double a11, double a12,
                                        double a21, double a22,
                                        double b1, double b2) const;

};


inline unsigned int
VoronoiCell::n_edges(void) const
{
  return(_edges.size());
}

inline const VoronoiCell::NodePair&
VoronoiCell::get_nodes(unsigned int e) const
{
  return(_edges[e]);
}

inline double
VoronoiCell::get_edge_length(unsigned int e) const
{
  return(_edge_lengths[e]);
}

inline const std::vector<double>&
VoronoiCell::get_edge_lengths(void) const
{
  return(_edge_lengths);
}

inline double
VoronoiCell::get_cell_face_area(unsigned int e) const
{
  return(_edge_areas[e]);
}

inline const std::vector<double>&
VoronoiCell::get_cell_face_areas(void) const
{
  return(_edge_areas);
}

inline double
VoronoiCell::get_volume(unsigned int n) const
{
  return(_volumes[n]);
}

inline const std::vector<double>&
VoronoiCell::get_volumes(void) const
{
  return(_volumes);
}


#endif // TC_VORONOICELL_H
