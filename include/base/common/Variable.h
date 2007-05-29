// $Id$

#ifndef _VARIABLE_H_
#define _VARIABLE_H_

#include "TypeDefs.h"

#include <map>
#include <string>


//! Interface for classes which contain variables
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
 * When the variable is used, it has to be referred to using the user defined
 * name \c variablename.
 *
 * A class derived from this interface can contain more than one variable
 * parameter. In this case, each such parameter is identified by an additional
 * id.
 */
class Variable
{

  private:

    //! The type for the list of variables
    typedef std::map<const std::string, std::pair<Variable*, ID> > VariableMap;


  public:

    //! Empty constructor
    Variable(void);


    //! Destructor
    virtual ~Variable(void);


    //! Get the variable object that has the variable \c var
    /*!
     * \param var the variable name
     * \return a pointer to the variable object or NULL, if \c var doesn't exist
     */
    static Variable* get_variable(const std::string& var);


    //! Check if \c var is a valid variable variable
    /*!
     * \param var the name of the variable
     * \return \c true, if \c var is defined, \c false otherwise
     */
    static bool is_variable(const std::string& var);


    //! Set the value of variable \c var
    /*!
     * \param var the variable name
     * \param value the value to set
     */
    static void set_variable_value(const std::string& var, double value);

    
    //! Get the value of variable \c var
    /*!
     * \param var the variable name
     * \param value the value to set
     */
    static double get_variable_value(const std::string& var);

    

  protected:

    //! Parse the string from the input file and register if needed
    /*!
     * A variable variable is defined in the input file using the notation
     * \c \#name(defaultvalue)
     *
     * If a string of this type is found, the variable object is registered
     * in the list and the default value is returned.
     *
     * \param s the string from the input file
     * \param defaultval the default value if not given in input file
     * \param id the class internal ID of the variable. This makes possible to have
     * more than one variable in a class
     * \return the default value as given in the inputfile
     */
    double check_and_register(const std::string& s, double defaultval, ID id = 0);


    //! Set the value of variable with ID \c id
    virtual void set_variable_value(double value, ID id = 0) = 0;


    //! Get the value of variable with ID \c id
    virtual double get_variable_value(ID id = 0) = 0;



  private:


    //! A list with all models that define a variable variable
    /*!
     * A variable variable is defined in the input file using the notation
     * \c @name(defaultvalue)
     */
    static VariableMap _variables;

};



//
// inline members
// 

inline 
Variable::Variable(void)
{
}


inline
Variable::~Variable(void)
{
}



#endif // _VARIABLE_H_
