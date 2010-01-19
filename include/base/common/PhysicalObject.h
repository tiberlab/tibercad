// $Id$

#ifndef _PHYSICALOBJECT_H_
#define _PHYSICALOBJECT_H_

#include "TiberModelObject.h"

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
     * Add a new \c PhysicalModelInterface object
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


    //! Get the material name
    const std::string& get_name(void) const;



  protected:

    //! a typedef for convenience
    typedef std::map<ID, PhysicalModel*> ModelMap;


    //! The constructor
    /*!
     * \param type the object type (BULK, BOUNDARY, EDGE, NODE)
     */
    PhysicalObject(ObjectType type, const ModelOptions& options);


    //! Set the objects name
    void set_name(const std::string& name);


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

    //! The objects name
    std::string _name;


    //! The database to be used
    Database* _database;


    //! The type of this object
    ObjectType _type;


    //! A flag to tell if the material is already initialized
    bool _is_initialized;


    //! The map containing all \c PhysicalModelInterface objects
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
const std::string&
PhysicalObject::get_name(void) const
{
  return  _name;
}


inline
void
PhysicalObject::set_name(const std::string& name)
{
  _name = name;
}


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

#endif // _PHYSICALOBJECT_H_
