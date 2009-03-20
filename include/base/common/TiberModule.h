// $Id$

#ifndef _TIBERMODULE_H_
#define  _TIBERMODULE_H_


#ifdef BUILD_TIBER_MODULES

//
// Provides macros needed to create a shared TiberCAD module
//

#ifdef CYGWIN
#  define DLLEXPORT __declspec(dllexport)
#  define DLLLOCAL
#else
#  ifdef GCC_HASVISIBILITY
#    define TBDLEXPORT __attribute__ ((visibility("default")))
#    define TBDLLOCAL __attribute__ ((visibility("hidden")))
#  else
#    define TBDLEXPORT
#    define TBDLLOCAL
#  endif
#endif

/*!
 * \def TIBER_MODULE(classname, simname)
 *
 * \brief Creates methods to create and destroy a simulation object
 * 
 * In each implementation derived from SimulationInterface, put
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
    TBDLEXPORT void destroy(SimulationInterface* p) { \
      delete p; \
    } \
    TBDLEXPORT classname* create(void) { \
      return new classname(); \
    } \
    TBDLLOCAL const char* _tiber_module_ = #simname; \
    TBDLLOCAL const char* library_name(void) { \
      return _tiber_module_; \
    } \
  }



/*!
 * \def TIBER_SUBMODEL(classname, modname)
 *
 * \brief Creates methods to create and destroy a simulation object
 * 
 * In each implementation derived from SimulationInterface, put
 * this macro somewhere in the source file to be able to compile
 * it as TiberCad module.
 *
 * \param name the name of the class that should be 'creatable'
 * \param modname the name for this module
 *
 * \c modname will be used to create the library name, and the model
 * will have to be referred to in the input file by \c modname
 */
#define TIBER_SUBMODEL(classname, modname) \
  extern "C" { \
    TBDLEXPORT void destroy(PhysicalModelInterface* p) { \
      delete p; \
    } \
    TBDLEXPORT classname* create(void) { \
      return new classname(); \
    } \
    TBDLLOCAL const char* _tiber_module_ ## modname = #modname; \
    TBDLLOCAL const char* library_name(void) { \
      return _tiber_module_ ## modname; \
    } \
  }



#else

#define TBDLEXPORT
#define TBDLLOCAL
#define TIBER_MODULE(classname, simname)
#define TIBER_SUBMODEL(classname, modname)

#endif // BUILD_TIBER_MODULES


#endif // _TIBERMODULE_H_
