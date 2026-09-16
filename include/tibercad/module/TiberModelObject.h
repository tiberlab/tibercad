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
 * \file TiberModelObject.h
 * \brief Public tiberCAD API header.
 *
 * \note This header is part of the public tiberCAD API.
 *       API version: 3.5
 */


#ifndef TC_TIBERMODELOBJECT_H
#define TC_TIBERMODELOBJECT_H

#include "tibercad/base/ModelOptions.h"
#include "tibercad/base/Initializer.h"
#include "tibercad/base/tiber_dll.h"
#include "tibercad/base/libMeshDefs.h"

// For debugging
#include "libmesh/reference_counted_object.h"


// stringify
#ifndef xstr
#define xstr(a) stringify(a)
#endif
#ifndef stringify
#define stringify(a) #a
#endif


#ifndef MODULENAME
#define MODULENAME
#endif


template <typename T> class InitializerBase;


//! The base class for all TiberCAD model classes
/*!
 * TiberCAD model classes are classes derived from:
 * \li PhysicalModel
 * \li BoundaryProperties
 * \li SimulationInterface
 */
class TiberModelObject
  : public libMesh::ReferenceCountedObject<TiberModelObject>
{

  public:

    //! Destructor
    virtual ~TiberModelObject(void);

    //! Set options for this model
    /*!
     * The options are stored internally and are accessible through
     * special methods.
     * Options have to be specified at creation time.
     */
    void set_options(const ModelOptions& options, bool add = false);

    //! Get the options for this model
    const ModelOptions& get_options(void) const;


    //! Get the options for this model
    ModelOptions& get_options(void);


    //! Set the name of a model
    /*!
     * If for some reason the automatically assigned name from 
     * the input file is usnuitable, this method can be used
     * to reassign a new name.
     */
    void set_name(const std::string& name);


    //! Get the user defined name of this model
    /*!
     * This is the name of the model/simulations as read from the input file.
     */
    const std::string& get_name(void) const;


    //! Destroy an object
    /*!
     * \param p the pointer to the object to destroy
     */
    static void destroy(TiberModelObject* p);


  protected:


    //! The creation method signature
    typedef TiberModelObject* (*create_t)(const ModelOptions&);


    //! The destruction method signature
    typedef void (*destroy_t)(TiberModelObject*);


    //! The constructor
    /*!
     * \param options the options to be assigned to this object
     */
    explicit TiberModelObject(const ModelOptions& options);


    //! Try to create an object from a dynamic link library
    /*!
     * \param name the name of the library <c>lib</c><it>name</it><c>.so</c>
     * \param options the options to pass to the object (which could be used to
     *   create different objects depending on the options)
     *
     * \return \c nullptr if it could not find the library
     *
     */
    template <typename T>
    static T* create_from_library(const std::string& name,
        const ModelOptions& options = ModelOptions());


    //! Create an object from a given creator function
    static TiberModelObject* create_from_function(create_t create, destroy_t destroy,
        const ModelOptions& options = ModelOptions());


    //! Create as a clone from a given object
    /*! \param handle a pointer to hold any object needed in the class implementation
     *   of the creation method.
     */
    template <typename T>
    static T* create_from_object(const T* other);


    //! Tells if a parameter has been specified in the input file
    /*!
     * \param name the parameter name
     * \param override do or do not override from other places (if
     *   implemented)
     */
    bool has_parameter(const std::string& name,
        bool override = true) const;


    //! Tells if an option has been specified in the input file
    /*!
     * \param name the option name
     * \param override do or do not override from other places (if
     *   implemented)
     */
    bool has_option(const std::string& name,
        bool override = true) const;


    //! Get a parameter
    /*!
     * Parameters are basically the same as options (cf. get_option())
     * , but they can be variables.
     *
     * \param name the parameter name
     * \param variable the variable where the value will be put
     * \param override do or do not override from other places (if
     *   implemented)
     *
     * \par Example:
     *  Assume you have a class \c MyClass inherited from TiberModelObject
     *  which has a member \c _var that you want to be able to use as
     *  variable and which needs initialization each time you set its value.
     *  \code
     *  class MyClass : public TiberModelObject
     *  {
     *    public:
     *      void init(void);
     *
     *      void prepare(void);
     *
     *    private:
     *      int _var;
     *  };
     *
     *  void MyClass::prepare(void)
     *  {
     *    get_parameter("myParameter", _var, true, initializer(&MyClass::init));
     *  }
     *  \endcode
     */
    template <typename T>
    void get_parameter(const std::string& name, T& variable,
        bool override = true, InitializerBase<T>* initfunc = nullptr);


    //! Get a parameter which is a vector of values (of the same type)
    /*
     * \param name the name of the option
     * \param vec the vector, where the values will be stored. \c vec can
     * contain default values, but it's size will be changed according to
     * the vector found in the options.
     */
    template <typename T>
    void get_parameter(const std::string& name, std::vector<T>& vec,
        bool override = true);


    //! Get an option
    /*!
     * Options cannot be variables.
     *
     * \param name the parameter name
     * \param variable the variable where the value will be put
     * \param override do or do not override from other places (if
     *   implemented)
     */
    template <typename T>
    T get_option(const std::string& name, T default_value, bool override = true) const;


    //! Get a string option with default given as const char
    std::string get_option(const std::string& name,
        const char* default_value,
        bool override = true) const;


    //! Get an option which is a vector of values
    template <typename T>
    void get_option(const std::string& name, std::vector<T>& vec, bool override = true) const;


    //! Get an option which is a set of values
    template <typename T>
    void get_option(const std::string& name, std::set<T>& vec, bool override = true) const;


    //! Get an option which is a vector of values
    void get_option(const std::string& name, libMesh::RealVectorValue& vec,
        bool override = true) const;


    //! Get an option which is a vector of values
    void get_option(const std::string& name, libMesh::RealTensor& vec, bool override = true) const;


    //! Get an option as a Point
    void get_option(const std::string& name, Point& point,
        bool override = true) const;


    //! To override a parameter string from a strange location
    /*!
     * Normally parameters/options are read only from the local
     * ModelOptions object. Implement this method if you want to
     * override from other sources.
     *
     * \param name the name of the parameter
     * \param s the currently assigned string
     */
    virtual void override_parameter_string(const std::string& name,
        std::string& s) const;


    //! Create an initializer functor from a member function
    template <class C, typename T = double>
    InitializerBase<T>* initializer(void (C::*func)(void));


    //! Create an initializer functor from a member function
    template <class C, typename T>
    InitializerBase<T>* initializer(void (C::*func)(T&));



  private:


    //! The type for library handles
    typedef void* libhandle_t;


    //! Don't allow default constructor
    TiberModelObject(void) = delete;


    //! Don't allow copy construction
    TiberModelObject(const TiberModelObject& other) = delete;


    //! Don't allow assignment operator
    TiberModelObject& operator=(const TiberModelObject& rhs) = delete;


    //! Try to create an object from a dynamic link library
    /*!
     * \param name the name of the library <c>lib</c><it>name</it><c>.so</c>
     * \param options the options to pass to the object (which could be used to
     * \return \c NULL if library cannot be opened
     */
    static TiberModelObject* _create_from_library(const std::string& name,
        const ModelOptions& options) TC_DLLOCAL;


    //! Try to create as a clone from an existing object
    /*! \param handle a pointer to hold any object needed in the class implementation
     *   of the creation method.
     */
    static TiberModelObject* _create_from_object(const TiberModelObject* other) TC_DLLOCAL;


    //! The options for this model as read from the input file
    ModelOptions _options;


    //! The library handle for this object
    libhandle_t _libhandle;


    //! The creation method for this object
    create_t _create;


    //! The destruction method for this object
    destroy_t _destroy;


    //! A user defined name for this model
    /*!
     * The name is assigned from the ModelOptions.
     */
    std::string _name;

};



//
// inline members
//



inline
void
TiberModelObject::set_options(const ModelOptions& options, bool add)
{
  if (add)
    _options += options;
  else
    _options = options;
}


inline
bool
TiberModelObject::has_option(const std::string& name, bool override) const
{
  return has_parameter(name, override);
}



inline
const ModelOptions&
TiberModelObject::get_options(void) const
{
  return _options;
}


inline
ModelOptions&
TiberModelObject::get_options(void)
{
  return _options;
}




inline
const std::string&
TiberModelObject::get_name(void) const
{
  return _name;
}



inline
void
TiberModelObject::set_name(const std::string& name)
{
  _name = name;
}



inline
void
TiberModelObject::override_parameter_string(const std::string&,
        std::string&) const
{
}





template <class C, typename T>
inline
InitializerBase<T>*
TiberModelObject::initializer(void (C::*func)(void))
{
  return new Initializer<C, T>(static_cast<C*>(this), func);
}


template <class C, typename T>
inline
InitializerBase<T>*
TiberModelObject::initializer(void (C::*func)(T&))
{
  return new Initializer<C, T>(static_cast<C*>(this), func);
}


template <typename T>
inline
T*
TiberModelObject::create_from_library(const std::string& name,
    const ModelOptions& options)
{
#ifdef TC_BUILD_TIBER_MODULES
  return dynamic_cast<T*>(_create_from_library(name, options));
#else
  return nullptr;
#endif
}

template <typename T>
inline
T*
TiberModelObject::create_from_object(const T* other)
{
#ifdef TC_BUILD_TIBER_MODULES
  return dynamic_cast<T*>(_create_from_object(other));
#else
  return nullptr;
#endif
}

#endif /* _TIBERMODELOBJECT_H_ */
