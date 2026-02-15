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
 * \file tiber_dll.h
 * \brief Public tiberCAD API header.
 *
 * \note This header is part of the public tiberCAD API.
 *       API version: 3.5
 */


#ifndef _TIBER_DLL_H_
#define _TIBER_DLL_H_

#include "tibercad/base/tiber_config.h"

#if defined(_WIN32)
# define TBDLEXPORT __declspec(dllexport)
# define TBDLLOCAL __declspec(dllimport)
#else
# ifdef TC_HAVE_FUNC_ATTRIBUTE_VISIBILITY
#   define TBDLEXPORT __attribute__ ((visibility("default")))
#   define TBDLLOCAL __attribute__ ((visibility("hidden")))
# else
#   define TBDLEXPORT
#   define TBDLLOCAL
# endif
#endif


// these are the symbol names for the dll entry points
#define TBCREATEFUNC tc_model_create
#define TBDESTROYFUNC tc_model_destroy
#define TBCREATEFUNCSYM "tc_model_create"
#define TBDESTROYFUNCSYM "tc_model_destroy"

#endif // _TIBER_DLL_H_
