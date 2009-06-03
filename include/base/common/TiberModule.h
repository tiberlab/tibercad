// $Id$

#ifndef _TIBERMODULE_H_
#define _TIBERMODULE_H_



//
// Provides macros needed to create a shared TiberCAD module
//

#ifdef BUILD_TIBER_MODULES
# ifdef CYGWIN
#  define DLLEXPORT __declspec(dllexport)
#  define DLLLOCAL
# else
#  ifdef GCC_HASVISIBILITY
#    define TBDLEXPORT __attribute__ ((visibility("default")))
#    define TBDLLOCAL __attribute__ ((visibility("hidden")))
#  else
#    define TBDLEXPORT
#    define TBDLLOCAL
#  endif
# endif

/*!
 * \def TIBER_MODULE(classname, simname)
 *
 * \brief Creates methods to create and destroy a simulation object
 *
 * In each implementation derived from TiberModelObject, put
 * this macro somewhere in the source file to be able to compile
 * it as TiberCad module.
 *
 * \param name the name of the class that should be 'creatable'
 * \param simname the name for this module
 *
 * \c simname will be used to create the library name, and the model
 * will have to be referred to in the input file by \c simname
 */
#define TIBER_MODULE(classname, simname) \
  extern "C" { \
    TBDLEXPORT void destroy(TiberModelObject* p) { \
      delete p; \
    } \
    TBDLEXPORT classname* create(void) { \
      return new classname(); \
    } \
  } \


#else

# define TBDLEXPORT
# define TBDLLOCAL
# define TIBER_MODULE(classname, simname)

#endif // BUILD_TIBER_MODULES


#endif // _TIBERMODULE_H_
