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


#ifndef _TIBERMODULE_H_
#define _TIBERMODULE_H_

#include "tibercad/base/tiber_dll.h"


//
// Provides macros needed to create a shared TiberCAD module
//



#ifndef MODULE_NAME
#define MODULE_NAME
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
#ifdef CREATABLE

#ifndef xstr
#define xstr(a) stringify(a)
#endif
#ifndef stringify
#define stringify(a) #a
#endif

extern "C" {
  TBDLEXPORT void
  TBDESTROYFUNC(TiberModelObject* p) {
    delete p;
  }

  TBDLEXPORT TiberModelObject*
  TBCREATEFUNC(const ModelOptions& options, const void* handle) {
    TiberModelObject* obj = NULL;
#ifdef CREATORCODE
#include xstr(CREATORCODE)
#else
    obj = CREATABLE::create(options);
    //obj = new CREATABLE(options); maybe should change to this
#endif
    static_cast<const void*>(handle);
    return obj;
  }
}
#endif




#endif // _TIBERMODULE_H_
