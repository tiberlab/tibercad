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
 * \file SpaceTransformation.h
 * \brief Public tiberCAD API header.
 *
 * \note This header is part of the public tiberCAD API.
 *       API version: 3.5
 */



#ifndef TC_SPACETRANSFORMATION_H
#define TC_SPACETRANSFORMATION_H

#include "tibercad/base/libMeshDefs.h"

#include "libmesh/vector_value.h"
#include "libmesh/tensor_value.h"

#include <vector>


/*!
 * \brief Implements space transformations
 *
 * This class implements all necessary space transformations.
 * In particular, rotations around arbitrary axes, mirroring at
 * arbitrary planes, and all elements of the point groups. The
 * latter are defined in canonical space coordinates x, y, z.
 *
 * Specific information has been taken from
 * http://lampx.tugraz.at/~hadley/ss2/crystalphysics/crystalclasses/crystalclasses.html
 *
 */
class SpaceTransformation
{

  public:

    /*! \brief Rotate around specified axis 
     *
     * Rotate a Point around a specified arbitrary axis, by
     * given angle.
     *
     * \param axis the rotation axis
     * \param angle the rotation angle in radians
     * \param point the point to be rotated
     */
    static void rotate(const libMesh::RealVectorValue& axis,
                       const double angle,
                       libMesh::Point& point);


    /*! \brief Mirror at specified plane
     *
     * Mirror a point at a specified plane
     *
     * \param normal the plane normal
     * \param origin a point on the plane
     * \param point the point to be mirrored
     */
    static void mirror(const libMesh::RealVectorValue& axis,
                       const libMesh::Point& origin,
                       libMesh::Point& point);


    /*! \brief Return all equivalent points according to given point group
     *
     * This method first (if necessary) calculates all symmetry operations
     * from the generators, and then constructs all equivalent points
     * for the provided input point. The resulting vector does not contain
     * duplicate points. Operations are done in reference coordinates.
     *
     * \param symmetry the space group (in any notation)
     * \param point original point
     * \param star the vector containing all equivalent points
     */
    static void create_star(const std::string& symmetry,
                            const libMesh::Point& point,
                            std::vector<libMesh::Point>& star);


  private:

    /*! \brief Constructor
     *
     * All operations are using static members
     */
    SpaceTransformation(void);

    /*! \brief Generate all transformation matrices for given symmetry
     *
     * \param symmetry Symmetry as International Symbol
     */
    static void generate_transformations(const std::string& symmetry);


};


#endif // TC_SPACETRANSFORMATION_H
