// $Id$

#ifndef _TIBERMODELOBJECT_H_
#define _TIBERMODELOBJECT_H_

#include "ModelOptions.h"
#include "Initializer.h"
#include "TypeDefs.h"

// For debugging
#include "reference_counted_object.h"


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

    //! Get the options for this contact model
    const ModelOptions& get_options(void) const;


    //! Get the options for this contact model
    ModelOptions& get_options(void);



  protected:

    //! Default constructor
    TiberModelObject(void) {};

    //! Copy constructor
    TiberModelObject(const TiberModelObject& other);


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
     *  /endcode
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

    //! Don't allow assignment operator
    TiberModelObject& operator=(const TiberModelObject& rhs);


    //! The options for this model as read from the input file
    ModelOptions _options;

};



//
// inline members
//

inline
TiberModelObject::TiberModelObject(const TiberModelObject& other)
  : _options(other._options)
{
}


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
void
TiberModelObject::override_parameter_string(const std::string& name,
        std::string& s) const
{
  ignore_unused_variable(name);
  ignore_unused_variable(s);
}



template <typename T>
inline
InitializerBase*
TiberModelObject::initializer(void (T::*func)(void))
{
  return new Initializer<T>(*static_cast<T*>(this), func);
}


#endif /* _TIBERMODELOBJECT_H_ */
