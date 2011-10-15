// $Id$

#ifndef _TIBERMODULE_H_
#define _TIBERMODULE_H_

#include "tiber_dll.h"


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
#endif
    return obj;
  }
}
#endif




#endif // _TIBERMODULE_H_
