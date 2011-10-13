// $Id$

#ifndef _TIBERMODULE_H_
#define _TIBERMODULE_H_

#include "tiber_dll.h"


//
// Provides macros needed to create a shared TiberCAD module
//

#define TBCREATEFUNC __create
#define TBDESTROYFUNC __destroy
#define TBCREATEFUNCSYM "__create"
#define TBDESTROYFUNCSYM "__destroy"


#ifndef MODULE_NAME
#define MODULE_NAME
#endif

/*!
 * \def TIBER_MODULE
 *
 * \brief Creates methods to create and destroy a simulation object
 *
 * In each implementation derived from TiberModelObject, put
 * this macro somewhere in the source file to be able to compile
 * it as TiberCad module.
 *
 */
#ifdef CREATABLE
extern "C" {
  TBDLEXPORT void
  TBDESTROYFUNC(TiberModelObject* p) {
    delete p;
  }

  TBDLEXPORT TiberModelObject*
  TBCREATEFUNC(const ModelOptions& options, const void* handle) {
    TiberModelObject* obj = NULL;
#ifdef CREATORCODE
#include CREATORCODE
#else
    obj = CREATABLE::create(options);
#endif
    return obj;
  }
}
#endif




#endif // _TIBERMODULE_H_
