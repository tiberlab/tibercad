// $Id$

#ifndef _VARIABLE_H_
#define _VARIABLE_H_


#include <map>
#include <string>

class TiberModelObject;
class InitializerBase;

//! Interface for variables
/*!
 * A variable is meant to be a parameter that can be changed at runtime,
 * eg. for doing a parameter sweep.
 *
 * In the input file, a variable has to be defined as follows:
 * \verbatim parameter = @variablename[defaultvalue] \endverbatim
 * where \c parameter is the identifier of the model parameter,
 * \c variablename is the user defined name for the variable and \c defaultvalue
 * is the default value to be assigned to this variable.
 *
 * A variable can be assigned to different parameteres in different
 * models. When assigning the value, it will be set in all objects that
 * use this variable.
 *
 * When the variable is used, it has to be referred to using the user defined
 * name \c variablename.
 *
 */
class Variable
{


  public:

    //! Destructor
    virtual ~Variable(void) { };


    //! Check if \c var is a valid variable variable
    /*!
     * \param var the name of the variable
     * \return \c true, if \c var is defined, \c false otherwise
     */
    static bool is_variable(const std::string& var);


    //! Get the name of this variable
    const std::string& get_name(void) const;


    //! Set the value of variable \c var
    /*!
     * \param var the variable name
     * \param value the value to set
     */
    template <typename T>
    static void set_variable_value(const std::string& var, const T& value);


    //! Get the value of variable \c var
    /*!
     * \param var the variable name
     * \return the value
     *
     * \note Makes no check about existence of the variable. Do this
     * first using \c is_variable()
     */
    template <typename T>
    static T get_variable_value(const std::string& var);


    //! Create a variable
    /*!
     * A variable is defined in the input file using the notation
     * \c @name(defaultvalue)
     *
     * If the variable already exists and has the same data type, the
     * given C++ variable is additionally assigned to it.
     *
     * \param s a string as read from the input file
     * \param variable the C++ variable
     * \param ct the container (class) holding \c variable
     *
     * \return a pointer to the
     * \c variable will contain the value read from the input file
     * if provided, else it will not change its value.
     */
    template <typename T>
    static void check_and_register(const std::string& s,
        T& variable, const TiberModelObject* ct = NULL,
        InitializerBase* initfunc = NULL);


    //! Unregister a model
    static void unregister(const TiberModelObject* ct);



  protected:

    //! Default constructor
    explicit Variable(const std::string& name);


    //! Do the unregistering of a model
    virtual void do_unregister(const TiberModelObject* ct) = 0;



  private:

    //! The type for the list of variables
    typedef std::map<const std::string, Variable*> VariableMap;


    //! The variable name
    const std::string _name;


    //! A list with all models that define a variable variable
    /*!
     * A variable variable is defined in the input file using the notation
     * \c @name(defaultvalue)
     *
     */
    static VariableMap _variables;

};



inline
const std::string&
Variable::get_name(void) const
{
  return _name;
}


#endif // _VARIABLE_H_
