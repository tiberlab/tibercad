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

class DofMap;
class MeshRegionInfo;
class BoundaryRegions;


class QuantumContact
{
  public:

    // Compute normal vector to a boundary region
    Point get_normal( unsigned int& count);

    // Extend mesh in normal vector direction
    void extend_mesh(void);

    // Project of any point in the extended mesh into a boundary point
    std::pair<const Elem*, Point> project_on_boundary (const Elem*,const Point&);

    //Calculate 3X3 matrix's determinant
    double Deter (const Point& P1, const Point& P2, const Point& P3);


    ID get_id(void) const;

    std::string get_name(void) const;

    const std::set<ID>& get_bd_ids(void) const;


    //! Destructor
    /*!
     * We do not declare it virtual here, as we will not allow
     * to derive from this class anyway.
     */
    ~QuantumContact(void);

    //! We need a public static creator function
    static QuantumContact* create(void);


    void init(const ID id,
              const std::string& name,
              Device* device,
              BoundaryRegions* bd_regions,
              const std::vector<ID>& rg_ids,
              const std::vector<ID>& bd_ids,
              double length);

    void plot(const std::string& name);

    void activate_elements(void);
    void inactivate_elements(void);

    //! this function must be called explicitly before use of qc
    //! since apparently the correct neighbor map gets broken
    void set_neighbor_map(void);

    static std::vector<ID> set2vec(const std::set<ID>& set);
    static std::set<ID> vec2set(const std::vector<ID>& vec);

  private:

    QuantumContact(void) TBDLLOCAL;

    void write_neighbors(void) const;

    Device* _device;

    MeshBase* _mesh;

    BoundaryRegions* _bd_regions;

    //! maps elements of the original mesh to neighboring elements in the quantum contact
    std::map<const Elem*, Elem*> _elemmap;

    //! maps elements of the quantum contact to side elements
    std::map<const Elem*,const ElementSide*> _elemsidemap;

    Point _normal;

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
std::string QuantumContact::get_name(void) const
{
  return _name;
}

#endif /* QUANTUMCONTACT_H_ */
