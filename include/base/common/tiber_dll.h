// $Id$

#ifndef _TIBER_DLL_H_
#define _TIBER_DLL_H_

#include "tiber_config.h"

#if defined(_WIN32)
#  define TBDLEXPORT __declspec(dllexport)
#  define TBDLLOCAL __declspec(dllimport)
# define TBDLEXPORT
# define TBDLLOCAL
#else
# ifdef GCC_HASVISIBILITY
#   define TBDLEXPORT __attribute__ ((visibility("default")))
#   define TBDLLOCAL __attribute__ ((visibility("hidden")))
# else
#   define TBDLEXPORT
#   define TBDLLOCAL
# endif
#endif


// these are the symbol names for the dll entry points
#define TBCREATEFUNC __create
#define TBDESTROYFUNC __destroy
#define TBCREATEFUNCSYM "__create"
#define TBDESTROYFUNCSYM "__destroy"

#endif // _TIBER_DLL_H_
