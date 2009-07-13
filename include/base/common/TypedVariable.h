// $Id$

#ifndef _TYPEDVARIABLE_H_
#define _TYPEDVARIABLE_H_

#include "Variable.h"
#include "InitializerBase.h"

#include <sstream>


//! The real variable class which is aware of the variable type
/*!
 * This class contains a mapping between C++ variable pointer
 * (usually a member variable) and the corresponding class which
 * has to be derived from TiberModelObject. It is possible to provide
 * the \c NULL pointer for the latter but <em> use this feature with
 * extreme care </em>, because if the C++ variable pointers gets
 * invalid it will \em not be unregistered.
 * You can additionally provide an initializer function which should
 * be called after every change of variable value. This method has to be a
 * member of the TiberModelObject class.
 */
template <typename T>
class TypedVariable : public Variable
{

  public:

    //! Constructor
    TypedVariable(const std::string& name)
      : Variable(name) { };

    //! Set the value
    void set_value(const T& value);

    //! Get the value
    const T& get_value(void) const { return _value; };

    //! Register a variable
    void register_variable(T& variable, const TiberModelObject* ct,
        InitializerBase* initfun);

    //! Get the value in string representation
    virtual std::string get_value_string(void) const;


  protected:

    /*! \copydoc Variable::do_unregister() */
    virtual void do_unregister(const TiberModelObject* ct);



  private:

    typedef std::pair<const TiberModelObject*, InitializerBase*> MapElem;
    typedef std::map<T*, MapElem> VarMap;

    //! The map containing all variables
    /*!
     * The void pointer points to the class containing the variable
     */
    VarMap _variables;


    //! The current value
    T _value;

};


template <typename T>
void
TypedVariable<T>::register_variable(T& variable, const TiberModelObject* ct,
    InitializerBase* initfunc)
{
  if (_variables.find(&variable) == _variables.end())
  {
    _variables[&variable] = MapElem(ct, initfunc);
    _value = variable;
  }
}


template <typename T>
void
TypedVariable<T>::set_value(const T& value)
{
  _value = value;
  typename VarMap::iterator it(_variables.begin());
  for ( ; it != _variables.end(); ++it)
  {
    *(it->first) = value;
    // call the initializer function
    if ((it->second).second != NULL)
      (*(it->second).second)();
  }
}



template <typename T>
void
TypedVariable<T>::do_unregister(const TiberModelObject* ct)
{
  typename VarMap::iterator it(_variables.begin());
  for ( ; it != _variables.end(); ++it)
  {
    if ((it->second).first == ct)
    {
      delete (it->second).second;
      _variables.erase(it);
    }
  }
}



template <typename T>
std::string
TypedVariable<T>::get_value_string(void) const
{
  std::ostringstream os;
  os << _value;

  return os.str();
}



#endif /* _TYPEDVARIABLE_H_ */
