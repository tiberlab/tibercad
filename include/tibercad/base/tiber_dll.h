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


#ifndef TC_TIBER_DLL_H
#define TC_TIBER_DLL_H

#include "tibercad/base/tiber_config.h"

#if defined(_WIN32)
# define TC_DLEXPORT __declspec(dllexport)
# define TC_DLLOCAL __declspec(dllimport)
#else
# ifdef TC_HAVE_FUNC_ATTRIBUTE_VISIBILITY
#   define TC_DLEXPORT __attribute__ ((visibility("default")))
#   define TC_DLLOCAL __attribute__ ((visibility("hidden")))
# else
#   define TC_DLEXPORT
#   define TC_DLLOCAL
# endif
#endif


// these are the symbol names for the dll entry points
#define TC_CREATEFUNC tc_model_create
#define TC_DESTROYFUNC tc_model_destroy
#define TC_CREATEFUNCSYM "tc_model_create"
#define TC_DESTROYFUNCSYM "tc_model_destroy"

#endif // TC_TIBER_DLL_H
