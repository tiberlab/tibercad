// $Id$

#ifndef _TIBER_DLL_H_
#define _TIBER_DLL_H_

#include "tiber_config.h"

#ifdef BUILD_TIBER_MODULES
# if defined(__CYGWIN__) || defined(__MINGW32__)
//#  define TBDLEXPORT __declspec(dllexport)
//#  define TBDLLOCAL __declspec(dllimport)
#  define TBDLEXPORT 
#  define TBDLLOCAL 
# else
#  ifdef GCC_HASVISIBILITY
#    define TBDLEXPORT __attribute__ ((visibility("default")))
#    define TBDLLOCAL __attribute__ ((visibility("hidden")))
#  else
#    define TBDLEXPORT
#    define TBDLLOCAL
#  endif
# endif
#else
# define TBDLEXPORT
# define TBDLLOCAL
#endif // BUILD_TIBER_MODULES


#endif // _TIBER_DLL_H_
