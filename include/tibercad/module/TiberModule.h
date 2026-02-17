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
 * \file TiberModule.h
 * \brief Public tiberCAD API header.
 *
 * \note This header is part of the public tiberCAD API.
 *       API version: 3.5
 */


#ifndef TC_TIBERMODULE_H
#define TC_TIBERMODULE_H

#include "tibercad/base/tiber_dll.h"


//
// Provides macros needed to create a shared TiberCAD module
//



#ifndef TC_MODULE_NAME
#define TC_MODULE_NAME
#endif

/*!
 * \brief Creates methods to create and destroy a simulation object
 *
 * In each implementation derived from TiberModelObject, you have
 * to include this header in the source file to be able to compile
 * it as TiberCad module.
 * For each module it may be included only once!
 *
 */
#ifdef TC_CREATABLE

template <class T>
class TC_ModelEnabler : public T
{
public:
    template <class... Args>
    TC_ModelEnabler(Args&&... args)
        : T(std::forward<Args>(args)...)
    {}
};

#ifndef TC_STRINGIFY
#define TC_STRINGIFY(a) TC_STRINGIFY1(a)
#define TC_STRINGIFY1(a) #a
#endif

extern "C" {
  TBDLEXPORT void
  TC_DESTROYFUNC(TiberModelObject* p) {
    delete p;
  }

  TBDLEXPORT TiberModelObject*
  TC_CREATEFUNC(const ModelOptions& options) {
    TiberModelObject* obj = NULL;
#ifdef TC_FACTORYCODE
#include TC_STRINGIFY(TC_FACTORYCODE)
#else
    obj = new TC_ModelEnabler<TC_CREATABLE>(options); 
#endif
    return obj;
  }
}

#endif // TC_CREATABLE




#endif // TC_TIBERMODULE_H
