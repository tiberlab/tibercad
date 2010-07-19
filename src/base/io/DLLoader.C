// $Id$

#include "tiber_config.h"
#include "DLLoader.h"
#include "Messages.h"
#include "TiberModule.h"
#include "RuntimeException.h"

#include <boost/filesystem/operations.hpp>

#include <dlfcn.h>
#include <sstream>

#ifdef CYGWIN
#define RTLD_NODELETE 0
#define RTLD_NOLOAD 0
#endif

//#ifdef DEBUG
//#define DLOPENFLAGS RTLD_NOW | RTLD_GLOBAL | RTLD_NODELETE
//#else
#define DLOPENFLAGS RTLD_LAZY | RTLD_GLOBAL | RTLD_NODELETE
//#endif


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
DLLoader::open_library(const string& name, DLLoader::LibraryInterface& iface)
{

  using namespace boost::filesystem;

  // we are not very optimistic
  bool success = false;

  // construct the library name
#ifdef CYGWIN
  string libfile = name + ".dll";
#else
  string libfile = name + ".so";
#endif

  bool file_exists = false;

  Messages::debug("Looking for library " + libfile + "... ");

  // we search for the library, as soon as we find it, we return
  // so it has the same behaviour as e.g. LD_LIBRARY_PATH
  list<string>::iterator it(_libpath.begin());
  const list<string>::iterator end(_libpath.end());
  for ( ; it != end; ++it)
  {
    path p(*it + "/" + libfile, native);
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
    //Messages::info("Trying to open " + libfile + "... ", false);

    const char* error_msg = 0;

    // check if it is already resident
    iface.handle = dlopen(libfile.c_str(), RTLD_NOLOAD | DLOPENFLAGS);

    if (iface.handle == NULL)
    {
      iface.handle = dlopen(libfile.c_str(), DLOPENFLAGS);

      // print the library name, but only the first time
      if ((iface.handle != NULL) && ((error_msg = dlerror()) == NULL))
        Messages::info("(using " + libfile + ")");
    }

    if (iface.handle != NULL)
    {
      iface.create_fnc = dlsym(iface.handle, TBCREATEFUNCSYM);
      iface.destroy_fnc = dlsym(iface.handle, TBDESTROYFUNCSYM);

      if ((iface.create_fnc == NULL) || (iface.destroy_fnc == NULL))
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
      os << "Error loading dynamic library " << libfile;
      throw RuntimeException(os.str());
    }
  }

}



  void
DLLoader::close_library(void* handle)
{
  if (handle != NULL)
    dlclose(handle);
}
