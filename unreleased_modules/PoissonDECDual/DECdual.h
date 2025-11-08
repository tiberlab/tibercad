#ifndef _DECDUAL_H_
#define _DECDUAL_H_

#include "OneChain2D.h"

#include <cstddef>
#include <vector>

namespace libMesh
{
  class Elem;
}

class PolyTope;

/*!
 * Representation of a dual element. The representation
 * is provided in terms of the center node and the star
 * of 1-chains. The duals of the latter form the boundary
 * of the dual element. Objects of this type rely on the
 * actual elements stored in a libMesh::Mesh object.
 * For now, the class is inteneded for 2D meshes only.
 * 
 * \brief Representation of a dual element.
 * \author Matthias Auf der Maur
 * \date 2025-03-20
 */
class DECdual
{

  public:

    /*!
     * Constructor
     * 
     * \param elem valid pointer to an element
     * \param node node index local to \c elem 
     */
    DECdual(const libMesh::Elem* elem, unsigned int node);

    DECdual(void) = delete;
    DECdual(DECdual&&) = default;
    DECdual(const DECdual&) = default;
    DECdual& operator=(DECdual&&) = default;
    DECdual& operator=(const DECdual&) = default;
    ~DECdual(void) = default;

    /*!
     * (Re)build the dual element from given geometry entities.
     * 
     * \param elem valid pointer to an element
     * \param node node index local to \c elem 
     * 
     * \brief Build a dual element
     */
    void build_dual_element(const libMesh::Elem* elem, unsigned int node);

    /*!
     * Get the Polytope associated to the dual element. It will be 
     * ordered counter-clock wise.
     * 
     * \param poly the PolyTope object to be filled
     */
    void get_polytope(PolyTope& poly) const;


  private:

    size_t _center_id;

    std::vector<OneChain2D> _star;

};


#endif // _DECDUAL_H_
