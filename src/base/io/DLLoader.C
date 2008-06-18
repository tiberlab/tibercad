// $Id$

#include "DLLoader.h"

#include <boost/filesystem/operations.hpp>

#include <dlfcn.h>
#ifdef DEBUG
#include <iostream>
#endif


using namespace std;


list<string>
DLLoader::_libpath;



void
DLLoader::set_library_path(const std::string& path)
{
  _libpath.clear();
  _libpath.push_front(path);
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




bool
DLLoader::open_library(const string& name, DLLoader::LibraryInterface& iface)
{

  using namespace boost::filesystem;

  // we are not very optimistic
  bool success = false;
  
  // construct the library name
  string libfile = "lib" + name + ".so";

  bool file_exists = false;

#ifdef DEBUG_
  cerr << "Looking for library " + libfile + "... ";
#endif

  list<string>::iterator it(_libpath.begin());
  const list<string>::iterator end(_libpath.end());
  for ( ; it != end; ++it)
    if (exists(*it + "/" + libfile))
    {
      libfile = *it + "/" + libfile;
      file_exists = true;
      break;
    }


#ifdef DEBUG_
  if (file_exists)
    cerr << "found." << endl;
  else
    cerr << "not found." << endl;
#endif

  if (file_exists)
  {
#ifdef DEBUG_
    cerr << "Trying to open " + libfile + "... ";
#endif

    // we will set it to false if something bad happens
    success = true;

    const char* error_msg = 0;

    //iface.handle = dlopen(libfile.c_str(), RTLD_NOW);
    iface.handle = dlopen(libfile.c_str(), RTLD_LAZY);
    if ((iface.handle != NULL) && ((error_msg = dlerror()) == NULL))
    {
      iface.create_fnc = dlsym(iface.handle, "create");
      iface.destroy_fnc = dlsym(iface.handle, "destroy");

      if ((iface.create_fnc == NULL)
          || (iface.destroy_fnc == NULL)
          || (dlerror() != NULL))
        success = false;
    }
    else
      success = false;

#ifdef DEBUG_
    if (success)
      cerr << "OK" << endl;
    else
      if (error_msg != 0)
        cerr << "failed: " << error_msg << endl;
      else
        cerr << "failed" << endl;
#endif
  }

  return success;
}



  void
DLLoader::close_library(void* handle)
{
  if (handle != NULL)
    dlclose(handle);
}
