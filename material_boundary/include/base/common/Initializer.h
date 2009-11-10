// $Id$

#ifndef _INITIALIZER_H_
#define _INITIALIZER_H_

#include "InitializerBase.h"


//! The real initializer class aware of type
template <typename T>
class Initializer : public InitializerBase
{

  public:

    typedef void (T::*InitFunc)(void);

    //! The constructor
    Initializer(T& obj, InitFunc func)
      : InitializerBase(obj, static_cast<CastFunc>(func)) {};

    //! The operator
    void operator()(void);

};



template <typename T>
inline
void
Initializer<T>::operator()(void)
{
  if (_func != NULL)
    (static_cast<T*>(_obj)->*static_cast<InitFunc>(_func))();
}


#endif /* _INITIALIZER_H_ */
