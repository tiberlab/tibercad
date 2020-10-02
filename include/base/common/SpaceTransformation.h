// $Id$


#ifndef _SPACETRANSFORMATION_H_
#define _SPACETRANSFORMATION_H_

#include "libMeshDefs.h"

#include "libmesh/vector_value.h"

#include <vector>


/*!
 * \brief Implements space transformations
 *
 * This class implements all necessary space transformations.
 * In particular, rotations around arbitrary axes, mirroring at
 * arbitrary planes, and all elements of the point groups
 */
class SpaceTransformation
{

  public:

    /*! \brief Constructor
     *
     * Most of operations are using static members
     */
    SpaceTransformation(void);


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
     * \param symmetry the space group
     * \param point original point
     * \param star the vector containing all equivalent points
     */
    static void create_star(const std::string& symmetry,
                            const libMesh::Point& point,
                            std::vector<libMesh::Point>& star);


  private:

};


#endif // _SPACETRANSFORMATION_H_
