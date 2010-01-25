// $Id$

#ifndef _MESHREGIONINFO_H_
#define _MESHREGIONINFO_H_

#include "TypeDefs.h"
#include "HashMap.h"

#include <string>
#include <iostream>


//! A class containing information on the region IDs
/*!
 * Mesh region IDs are assumed to be unique, and the mapping with names
 * one-to-one.
 */
class MeshRegionInfo
{

  public:

    //! Destructor
    ~MeshRegionInfo(void);


    //! Associate a name to a ID
    /*!
     * Nonexistent IDs will not be added!
     */
    void set_name(ID id, const std::string& name);


    //! Get the ID for a given name
    ID get_id(const std::string& name) const;


    //! Set a ID
    /*!
     * A ID that is not yet registered will be added with empty name
     */
    void add_id(ID id);


    //! Check if an id is present
    bool has_id(ID id) const;


    //! Get the ID for a given name
    const std::string& get_name(ID id) const;


    //! Cleanup the structure
    void clear(void);


  private:

    //! ID to name map
    /*!
     * This map contains all IDs present in the mesh, even those
     * that have not a name assigned.
     */
    typedef TiberCad::HashMap<ID, std::string>::Type IDToNameMap;

    //! name to ID map
    typedef TiberCad::HashMap<std::string, ID>::Type NameToIDMap;


    //! The names corresponding to the IDs
    IDToNameMap _ids_to_names;

    //! The IDs corresponding to the names
    NameToIDMap _names_to_ids;


    //! the empty string
    const std::string _empty;
};


inline
void
MeshRegionInfo::set_name(ID id, const std::string& name)
{
  IDToNameMap::const_iterator it(_ids_to_names.find(id));
  if (it != _ids_to_names.end())
  {
    _ids_to_names[id] = name;
    _names_to_ids[name] = id;
  }
}


inline
ID
MeshRegionInfo::get_id(const std::string& name) const
{
  ID id = INVALID_ID;
  NameToIDMap::const_iterator it(_names_to_ids.find(name));
  if (it != _names_to_ids.end())
    id = it->second;

  return id;
}


inline
const std::string&
MeshRegionInfo::get_name(ID id) const
{
  IDToNameMap::const_iterator it(_ids_to_names.find(id));
  if (it != _ids_to_names.end())
    return it->second;
  else
    return _empty;
}


inline
void
MeshRegionInfo::add_id(ID id)
{
  _ids_to_names[id];
}


inline
bool
MeshRegionInfo::has_id(ID id) const
{
  return _ids_to_names.count(id);
}


inline
MeshRegionInfo::~MeshRegionInfo(void)
{
  clear();
}

inline
void
MeshRegionInfo::clear(void)
{
  _ids_to_names.clear();
  _names_to_ids.clear();
}



#endif // _REGIONINFO_H_
