// $Id$

#ifndef _TIBERMODULE_H_
#define _TIBERMODULE_H_

#include "tiber_dll.h"

#ifdef BUILD_TIBER_MODULES

//
// Provides macros needed to create a shared TiberCAD module
//

#define TBCREATEFUNC __create
#define TBDESTROYFUNC __destroy
#define TBCREATEFUNCSYM "__create"
#define TBDESTROYFUNCSYM "__destroy"

/*!
 * \def TIBER_MODULE(classname, model [, type])
 *
 * \brief Creates methods to create and destroy a simulation object
 *
 * In each implementation derived from TiberModelObject, put
 * this macro somewhere in the source file to be able to compile
 * it as TiberCad module.
 *
 * \param classname the name of the class that should be 'createable'
 * \param model the model family (e.g. recombination, mobility etc.)
 * \param type the specific model name, if applicable (e.g. srh, auger)
 *
 * \c model and \c type will be used to create the library name as
 * model
 */
#define TIBER_MODULE(classname, model, type...) \
  extern "C" { \
    TBDLEXPORT void TBDESTROYFUNC(TiberModelObject* p) { \
      delete p; \
    } \
    TBDLEXPORT classname* TBCREATEFUNC(const ModelOptions& options) { \
      return classname::create(options); \
    } \
  } \


#else

# define TIBER_MODULE(classname, model, ...)

#endif // BUILD_TIBER_MODULES

#ifndef MODULE_NAME
#define MODULE_NAME
#endif

#endif // _TIBERMODULE_H_
