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
 * \file DLLoader.h
 * \brief Private tiberCAD header.
 *
 * \internal
 */


#ifndef TC_DLLOADER_H
#define TC_DLLOADER_H

#include "tibercad/base/tiber_dll.h"

#include <list>
#include <string>

//! A wrapper class which tries to load a certain library
/*!
 * The given library is looked for in well defined paths. If it is found,
 * it is opened and well defined symbols are looked for.
 */
class TC_DLLOCAL DLLoader
{

  public:

    //! A structure which defines the library interface
    struct LibraryInterface
    {
      //! The library handle
      void* handle = nullptr;

      //! The creation method
      void* create_fnc = nullptr;

      //! The destroy method
      void* destroy_fnc = nullptr;

      LibraryInterface(void) = default;
    };

    
    //! Opens a library and gets the library interface
    /*!
     * \param name the name of the library (\em not the filename)
     * \param iface the LibraryInterface where the interface information will
     * be stored
     */
    static void open_library(const std::string& name, LibraryInterface& iface);

    
    //! Closes a library for the given handle
    static void close_library(void* handle);


    //! Set the library search path
    static void set_library_path(const std::string& path);


    //! Prepend to library search path
    static void prepend_to_library_path(const std::string& path);


    //! Append to library search path
    static void append_to_library_path(const std::string& path);

    //! Get the library search paths, colon separated
    static void get_library_path(std::string& path);


  private:

    //! For static use only
    DLLoader(void);

    //! The library search path
    static std::list<std::string> _libpath;

};



#endif // TC_DLLOADER_H
