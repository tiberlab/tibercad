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
 * \file PhysicalObject.h
 * \brief Public tiberCAD API header.
 *
 * \note This header is part of the public tiberCAD API.
 *       API version: 3.5
 */


#ifndef TC_PHYSICALOBJECT_H
#define TC_PHYSICALOBJECT_H

#include "tibercad/module/TiberModelObject.h"
#include "tibercad/base/TypeDefs.h"

#include <map>


class Database;
class PhysicalModel;


//! The base class for physical entities as bulk, sides, edges and nodes
class PhysicalObject : public TiberModelObject
{

  public:

    //! The object type (bulk, surface, edge, node)
    enum ObjectType
    {
      INVALID = 0,  //!< invalid type
      BULK,         //!< bulk, N dimensional (3D, 2D, 1D)
      BOUNDARY,     //!< boundary, N-1 dimensional (3D, 2D, 1D)
      EDGE,         //!< edge, N-2 dimensional (3D, 2D)
      NODE          //!< node, 0 dimensional (3D only)
    };

    //! Destructor
    /*!
     * Destroys all associated models
     */
    virtual ~PhysicalObject(void);

    
    //! Set the database to be used
    /*!
     * The database could in principle be overriden from the options
     */
    void set_database(const Database& database);


    //! Initialize the object
    /*!
     * Read all needed data from the database.
     * It calls the \c read_database() of each \c PhysicalModel
     * object.
     */
    void init(void);


    //! Add new physical model
    /*!
     * Add a new \c PhysicalModel object
     *
     * \param model the model to add
     * \param simulator_id the id of the simulator this model is used for
     */
    void add_model(PhysicalModel* model, ID simulator_id);


    //! Get physical model for simulator with ID id
    /*!
     *
     * It will return the \c NULL pointer if the requested model
     * is not in the list.
     * \c id is a unique ID assigned to each simulator.
     *
     * \param id the identifier for simulator
     * \return a pointer to the model object
     */
    PhysicalModel* get_model(ID id) const;


    //! Get a reference to the database
    const Database& get_database(void) const;


    //! Get a writable reference to the database
    Database& get_database(void);


    //! Get the type of this object
    ObjectType get_type(void) const;



  protected:

    //! a typedef for convenience
    typedef std::map<ID, PhysicalModel*> ModelMap;


    //! The constructor
    /*!
     * \param type the object type (BULK, BOUNDARY, EDGE, NODE)
     *
     * \note The options have to contain an option \c dimension
     * with the spatial dimension of this object!
     */
    PhysicalObject(ObjectType type, const ModelOptions& options);


    //! The real init function
    /*!
     * This one gets called from init()
     */
    virtual void do_init(void);


    //! Get an iterator to the first model
    ModelMap::iterator models_begin(void);


    //! Get an iterator to the last model
    ModelMap::iterator models_end(void);



  private:


    //! The database to be used
    Database* _database;


    //! The type of this object
    ObjectType _type;


    //! A flag to tell if the material is already initialized
    bool _is_initialized;


    //! The map containing all \c PhysicalModel objects
    /*!
     * This map containes the physical model of any simulation type
     * requested.
     */
    ModelMap _models;


};


//
// inline members
//




inline
PhysicalModel*
PhysicalObject::get_model(ID id) const
{
  const ModelMap::const_iterator end(_models.end());
  ModelMap::const_iterator it(_models.find(id));
  if (it != end)
    return it->second;
  else
    return NULL;
}


inline
const Database&
PhysicalObject::get_database(void) const
{
  return *_database;
}


inline
Database&
PhysicalObject::get_database(void)
{
  return *_database;
}




inline
PhysicalObject::ModelMap::iterator
PhysicalObject::models_begin(void)
{
  return _models.begin();
}


inline
PhysicalObject::ModelMap::iterator
PhysicalObject::models_end(void)
{
  return _models.end();
}


inline
PhysicalObject::ObjectType
PhysicalObject::get_type(void) const
{
  return _type;
}

#endif // TC_PHYSICALOBJECT_H
