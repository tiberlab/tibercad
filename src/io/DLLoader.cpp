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
 * \file DLLoader.C
 * \brief Internal tiberCAD code.
 *
 * \internal
 */


#include "DLLoader.h"
#include "tibercad/base/tiber_config.h"
#include "tibercad/io/Messages.h"
#include "tibercad/base/tiber_dll.h"
#include "tibercad/base/RuntimeException.h"

#include <boost/filesystem/operations.hpp>

#include <dlfcn.h>
#include <sstream>

#if defined(__CYGWIN__) || defined(__MINGW32__)
#define RTLD_NODELETE 0
#define RTLD_NOLOAD 0
#endif

#ifdef DEBUG
#define TC_DLOPENFLAGS RTLD_NOW | RTLD_LOCAL
#else
#define TC_DLOPENFLAGS RTLD_LAZY | RTLD_LOCAL
#endif


using namespace std;


list<string>
DLLoader::_libpath;



void
DLLoader::set_library_path(const std::string& path)
{
  _libpath.clear();
  _libpath.push_back(path);
}



void
DLLoader::prepend_to_library_path(const std::string& path)
{
  _libpath.push_front(path);
}



void
DLLoader::append_to_library_path(const std::string& path)
{
  _libpath.push_back(path);
}


void
DLLoader::get_library_path(std::string& path)
{
  path = "";
  bool first = true;
  for (auto &&p : _libpath)
  {
    if (!first)
      path += ":";

    path += p;

    first = false;
  }
}




void
DLLoader::open_library(const string& name, DLLoader::LibraryInterface& iface)
{

  using namespace boost::filesystem;


  // construct the library name
#if defined(__CYGWIN__) || defined(__MINGW32__)
  string libfile = name + ".dll";
#elif defined(__APPLE__)
  string libfile = name + ".dylib";
#else
  string libfile = name + ".so";
#endif

  bool file_exists = false;


  // we search for the library, as soon as we find it, we return
  // so it has the same behaviour as e.g. LD_LIBRARY_PATH
  list<string>::iterator it(_libpath.begin());
  const list<string>::iterator end(_libpath.end());
  for ( ; it != end; ++it)
  {
    path p(*it + "/" + libfile);
    Messages::debug("Looking for library " + p.string() + "... ");
    if (exists(p))
    {
      libfile = *it + "/" + libfile;
      file_exists = true;
      break;
    }
  }

  if (file_exists)
    Messages::debug("found.");
  else
    Messages::debug("not found.");


  if (file_exists)
  {
    // check if it is already resident
    iface.handle = dlopen(libfile.c_str(), RTLD_NOLOAD | TC_DLOPENFLAGS);

    if (iface.handle == nullptr)
    {
      iface.handle = dlopen(libfile.c_str(), TC_DLOPENFLAGS);

      if (iface.handle != nullptr)
        Messages::info("(using " + libfile + ")");
    }


    if (iface.handle != nullptr)
    {
      iface.create_fnc = dlsym(iface.handle, TC_CREATEFUNCSYM);
      iface.destroy_fnc = dlsym(iface.handle, TC_DESTROYFUNCSYM);

      if ((iface.create_fnc == nullptr) || (iface.destroy_fnc == nullptr))
      {
        ostringstream os;
        os << "Cannot use dynamic library " << libfile
            << " (missing creation/destruction methods)";
        throw RuntimeException(os.str());
      }
    }
    else
    {
      ostringstream os;
      os << "Error loading dynamic library " << libfile
          << " (" << dlerror() << ")";
      throw RuntimeException(os.str());
    }
  }

}



void
DLLoader::close_library(void* handle)
{
  if (handle != nullptr)
    dlclose(handle);
}
