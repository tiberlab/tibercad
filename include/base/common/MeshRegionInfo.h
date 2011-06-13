// $Id$

#ifndef _MESHREGIONINFO_H_
#define _MESHREGIONINFO_H_

#include "TypeDefs.h"
#include "HashMap.h"
#include "IDSet.h"
#include "tiber_dll.h"

#include <string>
#include <sstream>
#include <iostream>


//! A class containing information on the region IDs
/*!
 * Mesh region IDs are assumed to be unique, and the mapping with names
 * one-to-one.
 */
class TBDLLOCAL MeshRegionInfo
{

  public:

    //! Constructor
    MeshRegionInfo(void) {}

    //! Destructor
    virtual ~MeshRegionInfo(void);


    //! Associate a name to an ID
    /*!
     * Nonexistent IDs will not be added!
     */
    void set_name(ID id, const std::string& name);


    //! Get the IDs for a given name
    const IDSet& get_ids(const std::string& name) const;


    //! Set an ID
    /*!
     * A ID that is not yet registered will be added with the ID
     * as default name.
     */
    void add_id(ID id);


    //! Check if id is present
    bool has_id(ID id) const;


    //! Get the total number of regions
    unsigned int n_subdomains(void) const;


    //! Cleanup the structure
    void clear(void);


    //! Get the next valid ID
    ID next_id(void) const;


    void print_info() const;


  protected:


    //! Get the name for a given id
    const std::string& get_name(ID id) const;


  private:

    //! ID to name map
    /*!
     * This map contains all IDs present in the mesh, even those
     * that have not a name assigned.
     */
    typedef HashMap<ID, std::string>::Type IDToNameMap;

    //! name to ID map
    /*!
     * Each ID can in principle have different names associated.
     */
    typedef HashMap<std::string, std::set<ID> >::Type NameToIDMap;



    //! The names corresponding to the IDs
    IDToNameMap _ids_to_names;

    //! The IDs corresponding to the names
    NameToIDMap _names_to_ids;


    //! the empty string
    const std::string _empty_string;

    //! an empty set
    const IDSet _empty_set;
};


inline
void
MeshRegionInfo::set_name(ID id, const std::string& name)
{
  IDToNameMap::const_iterator it(_ids_to_names.find(id));
  if (it != _ids_to_names.end())
  {
    _ids_to_names[id] = name;
    _names_to_ids[name].insert(id);
  }
}


inline
const IDSet&
MeshRegionInfo::get_ids(const std::string& name) const
{
  NameToIDMap::const_iterator it(_names_to_ids.find(name));
  if (it != _names_to_ids.end())
    return it->second;

  return _empty_set;
}


inline
const std::string&
MeshRegionInfo::get_name(ID id) const
{
  IDToNameMap::const_iterator it(_ids_to_names.find(id));
  if (it != _ids_to_names.end())
    return it->second;
  else
    return _empty_string;
}


inline
void
MeshRegionInfo::add_id(ID id)
{
  if (_ids_to_names[id].size() == 0)
  {
    std::ostringstream os;
    os << id;
    _ids_to_names[id] = os.str();
    _names_to_ids[os.str()].insert(id);
  }
}


inline
bool
MeshRegionInfo::has_id(ID id) const
{
  return _ids_to_names.count(id);
}


inline
unsigned int
MeshRegionInfo::n_subdomains(void) const
{
  return _ids_to_names.size();
}



#endif // _REGIONINFO_H_
