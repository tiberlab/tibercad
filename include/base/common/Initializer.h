// $Id$

#ifndef _INITIALIZER_H_
#define _INITIALIZER_H_

#include "InitializerBase.h"


//! The real initializer class aware of type
template <class Class, typename T>
class TBDLLOCAL Initializer : public InitializerBase<T>
{

  public:

    typedef void (Class::*VoidInitFunc)(void);
    typedef void (Class::*InitFunc)(T&);

    //! The constructor
    Initializer(Class* obj, VoidInitFunc func)
      : _obj(obj), _func(NULL), _voidfunc(func) {};

    //! The constructor
    Initializer(Class* obj, InitFunc func)
      : _obj(obj), _func(func), _voidfunc(NULL) {};

    //! The operator
    void operator()(T& val);


  private:

    Class* _obj;

    InitFunc _func;

    VoidInitFunc _voidfunc;
};



template <class Class, typename T>
inline
void
Initializer<Class, T>::operator()(T& val)
{
  if (_func != NULL)
    (_obj->*_func)(val);
  else if (_voidfunc != NULL)
    (_obj->*_voidfunc)();
}

// We provide a few modifier Functors

//! Inversion
class Invert : public InitializerBase<double>
{
  public:
    void operator()(double& val) { val = 1.0 / val;}
};



#endif /* _INITIALIZER_H_ */
