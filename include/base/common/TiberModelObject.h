// $Id$

#ifndef _TIBERMODELOBJECT_H_
#define _TIBERMODELOBJECT_H_

#include "ModelOptions.h"
#include "Initializer.h"
#include "tiber_dll.h"

// For debugging
#include "reference_counted_object.h"


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


class InitializerBase;


//! The base class for all TiberCAD model classes
/*!
 * TiberCAD model classes are classes derived from:
 * \li PhysicalModelInterface
 * \li BoundaryProperties
 * \li SimulationInterface
 */
class TiberModelObject
  : public ReferenceCountedObject<TiberModelObject>
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
    void set_options(const ModelOptions& options);

    //! Get the options for this model
    const ModelOptions& get_options(void) const;


    //! Get the options for this model
    ModelOptions& get_options(void);


    //! Set the name of a model
    /*!
     * Use with caution as it could break standard behaviour!
     */
    void set_name(const std::string& name);


    //! Get the user defined name of this model
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
    TiberModelObject(const ModelOptions& options);


    //! Try to create an object from a dynamic link library
    /*!
     * \param name the name of the library <c>lib</c><it>name</it><c>.so</c>
     * \param options the options to pass to the object (which could be used to
     * create different objects depending on the options)
     *
     * \return \c NULL if it could not find the library
     *
     */
    template <typename T>
    static T* create_from_library(const std::string& name,
        const ModelOptions& options = ModelOptions()) TBDLLOCAL;


    //! Create an object from a given creator function
    static TiberModelObject* create_from_function(create_t create, destroy_t destroy,
        const ModelOptions& options = ModelOptions()) TBDLLOCAL;


    //! Create a new model of the same type
    virtual TiberModelObject* create_new(void) const;


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
     *  Assume you have a class \c Pippo inherited from TiberModelObject
     *  which has a member \c _var that you want to be able to use as
     *  variable and which needs initialization each time you set its value.
     *  \code
     *  class Pippo : public TiberModelObject
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
     *  void Pippo::prepare(void)
     *  {
     *    get_parameter("myParameter", _var, true, initializer(&Pippo::init));
     *  }
     *  \endcode
     */
    template <typename T>
    void get_parameter(const std::string& name, T& variable,
        bool override = true, InitializerBase* initfunc = NULL);


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
    T get_option(const std::string& name, T default_value,
        bool override = true) const;


    //! Get a string option with default given as const char
    std::string get_option(const std::string& name,
        const char* default_value,
        bool override = true) const;


    //! Get an option which is a vector of values
    template <typename T>
    void get_option(const std::string& name, std::vector<T>& vec,
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
    template <typename T>
    InitializerBase* initializer(void (T::*func)(void));



  private:


    //! The type for library handles
    typedef void* libhandle_t;


    //! Don't allow default constructor
    TiberModelObject(void) TBDLLOCAL;


    //! Don't allow copy construction
    TiberModelObject(const TiberModelObject& other) TBDLLOCAL;


    //! Don't allow assignment operator
    TiberModelObject& operator=(const TiberModelObject& rhs) TBDLLOCAL;


    //! Try to create an object from a dynamic link library
    /*!
     * \param name the name of the library <c>lib</c><it>name</it><c>.so</c>
     * \param options the options to pass to the object (which could be used to
     * \return \c NULL if library cannot be opened
     */
    static TiberModelObject* _create_from_library(const std::string& name,
        const ModelOptions& options = ModelOptions()) TBDLLOCAL;



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
TiberModelObject::set_options(const ModelOptions& options)
{
  _options += options;
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



template <typename T>
inline
InitializerBase*
TiberModelObject::initializer(void (T::*func)(void))
{
  return new Initializer<T>(*static_cast<T*>(this), func);
}


template <typename T>
inline
T*
TiberModelObject::create_from_library(const std::string& name,
    const ModelOptions& options)
{
#ifdef BUILD_TIBER_MODULES
  return dynamic_cast<T*>(_create_from_library(name, options));
#else
  return NULL;
#endif
}

#endif /* _TIBERMODELOBJECT_H_ */
