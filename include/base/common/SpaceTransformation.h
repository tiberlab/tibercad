// $Id$


#ifndef _SPACETRANSFORMATION_H_
#define _SPACETRANSFORMATION_H_

#include "libMeshDefs.h"

#include "libmesh/vector_value.h"


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
     * \param \c axis the rotation axis
     * \param \c angle the rotation angle in radians
     * \param \c point the point to be rotated
     */
    static void rotate(const libMesh::RealVectorValue& axis,
                       const double angle,
                       libMesh::Point& point);


    /*! \brief Mirror at specified plane
     *
     * Mirror a point at a specified plane
     *
     * \param \c normal the plane normal
     * \param \c origin a point on the plane
     * \param \c point the point to be mirrored
     */
    static void mirror(const libMesh::RealVectorValue& axis,
                       const libMesh::Point& origin,
                       libMesh::Point& point);


  private:

};


#endif // _SPACETRANSFORMATION_H_
