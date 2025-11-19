// $Id$

#ifndef _DATACACHE_H
#define _DATACACHE_H

#include "tibercad/base/TypeDefs.h"
#include "tibercad/base/libMeshDefs.h"

#include <map>
#include <vector>
#include <unordered_map>
#include <set>
#include <deque>

#include "point.h"

/*!
 *
 * \brief A cache for solution values on mesh points
 *
 * This container serves to cache data requested from modules
 * in order to increase efficiency. This is useful whenever different models
 * access the same data from another module. Moreover, it allows to handle
 * circular dependencies.
 * In the current implementation, maximum cache size is not enforced.
 */
class DataCache
{

  public:

    //! Constructor
    /*!
     * The depth indicates how many elements are stored
     */
    DataCache(size_t size = 0);

    //! Destructor
    ~DataCache(void);

    //! Get data on a point
    /*!
     * \return false if data is not present
     * \input p a point in local coordinates
     *
     * Data will be returned in the input reference \c data
     */
    bool get_data(const libMesh::Elem* elem, const libMesh::Point& p,
        ID id, std::vector<double>& data) const;
    
    //! Get data on a series of points
    /*!
     * \return a set with the vector indices of the points found,
     *    referring to the input vector \c points
     * \input points a vector with points in local coordinates
     * \input data holds the data in the same layout as used
     *    in the \c get_solution() methods of \c SimulationInterface
     */
    std::set<unsigned int> get_data(const libMesh::Elem* elem,
        const std::vector<libMesh::Point>& points,
        ID id, std::vector<double>& data);

    //! Add data to the cache
    /*!
     * \input p a point in local coordinates
     * \input id the ID of the solution variable
     * \input data the associated data
     */
    void put_data(const libMesh::Elem* elem, const libMesh::Point& p,
        ID id, const std::vector<double>& data);


    //! Add data for several points to cache
    void put_data(const libMesh::Elem* elem,
        const std::vector<libMesh::Point>& points,
        std::map<ID, std::vector<double>>& data);


    //! Erase whole cache
    void flush(void);

    //! Get the memory occupation in Bytes
    size_t get_memory_size(void) const;


  private:

    //! Hash function for point coordinates, relying on being local coordinates
    class PointHash
    {
      public:
        size_t operator()(const libMesh::Point& p) const
        {
          std::hash<size_t> h;
          return(h(p(0)+10*p(1)+100*p(2)));
        }
    };

    //! The ID to data map
    typedef std::map<ID, std::vector<double>> Data;

    //! The data map type
    typedef std::unordered_map<libMesh::Point, std::deque<Data>::iterator, PointHash> PointMap;

    //! The external (element) container type
    typedef std::unordered_map<const libMesh::Elem*, PointMap> ElemMap;

    //! The cached elements
    mutable ElemMap _elements;

    //! The actual data
    std::deque<Data> _data;

    //! The cache size
    size_t _size;

    //! The currently accessed element
    mutable ElemMap::iterator _current;

};


#endif // _DATACACHE_H
