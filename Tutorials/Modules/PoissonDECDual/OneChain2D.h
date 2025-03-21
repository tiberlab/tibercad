#ifndef _ONECHAIN2D_H_
#define _ONECHAIN2D_H_

namespace libMesh
{
  class Elem;
}

/*!
 * A representation of a 1-chain in 2D. Internally, it
 * is stored as a pair of an element counter and a side
 * number. This allows to obtain the 1-chain itself,
 * i.e. the element side (a 1D segment) and its dual
 * (the segment between the element barycenters).
 *
 * \brief A representation of a 1-chain in 2D
 * \author Matthias Auf der Maur
 * \date 2025-03-20
 */
class OneChain2D
{
  public:

    OneChain2D(void) = delete;
    OneChain2D(OneChain2D&&) = default;

    /*!
     * The constructor
     * \param elem the element pointer
     * \param side the side number
     * 
     * \brief The constructor to be used.
     */
    OneChain2D(const libMesh::Elem* elem, unsigned int side);

    //! Get the element pointer
    const libMesh::Elem* get_elem(void) const;

    //! Get the side index
    unsigned int get_side(void) const;

  private:

    //! The element pointer
    const libMesh::Elem *_elem = nullptr;

    //! The side number
    unsigned int _side = 0;
};


inline
const libMesh::Elem*
OneChain2D::get_elem(void) const
{
  return(_elem);
}

inline
unsigned int
OneChain2D::get_side(void) const
{
  return(_side);
}


#endif // _ONECHAIN2D_H_
