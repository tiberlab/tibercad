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
 * \file Initializer.h
 * \brief Public tiberCAD API header.
 *
 * \note This header is part of the public tiberCAD API.
 *       API version: 3.5
 */


#ifndef _INITIALIZER_H_
#define _INITIALIZER_H_

#include "tibercad/base/InitializerBase.h"


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
