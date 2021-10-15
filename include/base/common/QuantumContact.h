/*
 * QuantumContact.h
 *
 *  Created on: Nov 3, 2011
 *      Author: fpalomba
 */

#ifndef QUANTUMCONTACT_H_
#define QUANTUMCONTACT_H_

#include "SimulationInterface.h"
#include "ModelOptions.h"
#include "ElementSide.h"


#include <string>

namespace libMesh
{
  class DofMap;
}

class BoundaryRegions;

/*
 * A QuantumContact is extruded from a set of boundaries of a set regions (reg_ids)
 *
 * Problem: currently we assign a unique id to a QuantumContact and to all its elem.
 *          
 *          this is not good in case we extrude from different material regions.
 *
 * Extrusion length for each of the two layers can be provided from the input. If it
 * is not specified, it will be calculated internally based on the mean height of the
 * elements at the surface. Thei height is obtained by normal projection of the points
 * not on the surface.
 *
 */
class QuantumContact
{
  public:

    typedef std::map<const libMesh::Elem*, const ElementSide*> ContactElemMap;
    typedef std::map<const libMesh::Elem*, const ElementSide*>::iterator ContactElemIterator;

    // get normal 
    const Point& get_normal(void) const;
    
    // get area
    double get_area(void) const;

    // Extend mesh in normal vector direction
    void extend_mesh(void);

    // Project of any point in the extended mesh into a boundary point
    std::pair<const libMesh::Elem*, libMesh::Point> project_on_boundary(const libMesh::Elem*, const Point&);

    // Project of point in the extended mesh into a boundary point
    std::pair<const libMesh::Elem*, std::vector<libMesh::Point> > project_on_boundary(const libMesh::Elem*, 
                                                  const std::vector<libMesh::Point>& p);


    ID get_id(void) const;

    std::string get_name(void) const;

    const std::set<ID>& get_bd_ids(void) const;

    //! Get the internal region ids touching the contact
    const std::set<ID>& get_region_ids(void) const;

    //! Get the iterator to the first contact element
    ContactElemIterator contact_elements_begin(void);

    //! Get the past-the-end iterator for the contact elements
    ContactElemIterator contact_elements_end(void);

    //! Destructor
    /*!
     * We do not declare it virtual here, as we will not allow
     * to derive from this class anyway.
     */
    ~QuantumContact(void);

    //! We need a public static creator function
    static QuantumContact* create(void) TBDLLOCAL;


    /*!
     * \brief Initialize a new QuantumContact object
     *
     * \param id e new unique region ID
     * \param name a unique name (read from input file)
     * \param device the device it is associated with
     * \param bd_regions the boundary regions container
     * \param rg_ids the IDs of the regions touched
     * \param bd_ids the IDs of the boundaries it touches
     * \param length the extrusion length
     *
     * If length = 0, it will be calculated internally as
     * the mean radius of the elements touching the contact.
     */
    void init(const ID id,
              const std::string& name,
              Device* device,
              BoundaryRegions* bd_regions,
              const std::set<ID>& rg_ids,
              const std::set<ID>& bd_ids,
              double length = 0.0)  TBDLLOCAL;

    void plot(const std::string& name);

    //! this function must be called explicitly before use of qc
    //! since apparently the correct neighbor map gets broken
    void set_neighbor_map(void);

    static std::vector<ID> set2vec(const std::set<ID>& set);
    static std::set<ID> vec2set(const std::vector<ID>& vec);

  private:

    QuantumContact(void) TBDLLOCAL;

    //! Compute normal vector to a boundary region
    Point get_normal( double& area);

    void write_neighbors(void) const;

    //Calculate 3X3 matrix's determinant
    double Deter(const Point& P1, const Point& P2, const Point& P3);

    Device* _device;

    MeshBase* _mesh;

    BoundaryRegions* _bd_regions;

    //! maps elements of the original mesh to neighboring elements in the quantum contact
    std::map<const Elem*, Elem*> _elemmap;

    //! maps elements of the quantum contact to side elements
    std::map<const Elem*,const ElementSide*> _elemsidemap;

    Point _normal;

    double _area;

    ID _id;

    std::string _name;

    std::set<ID> _rg_ids;

    std::set<ID> _bd_ids;

    double _length;

};

inline
ID QuantumContact::get_id(void) const
{
  return _id;
}

inline
const std::set<ID>& QuantumContact::get_bd_ids(void) const
{
  return _bd_ids;
}

inline
const std::set<ID>& QuantumContact::get_region_ids(void) const
{
  return _rg_ids;
}

inline
std::string QuantumContact::get_name(void) const
{
  return _name;
}

inline
const Point& QuantumContact::get_normal(void) const
{
   return _normal;
}

inline
double QuantumContact::get_area(void) const
{
   return _area;
}

inline
QuantumContact::ContactElemIterator
QuantumContact::contact_elements_begin(void)
{
  return(_elemsidemap.begin());
}

inline
QuantumContact::ContactElemIterator
QuantumContact::contact_elements_end(void)
{
  return(_elemsidemap.end());
}

#endif /* QUANTUMCONTACT_H_ */
