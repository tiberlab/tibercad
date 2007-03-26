// $Id$

#include "DLLoader.h"


#include <dlfcn.h>
#ifdef DEBUG
#include <iostream>
#endif


using namespace std;


string
DLLoader::_libpath = ".";

bool
DLLoader::open_library(const string& name, DLLoader::LibraryInterface& iface)
{

  bool success = true;
  
  string libfile = _libpath + "/" + "lib" + name + ".so";

#ifdef DEBUG_
  cerr << "Looking for library " + libfile + "... ";
  cerr << "found." << endl;
#endif

#ifdef DEBUG_
  cerr << "Trying to open " + libfile + "... ";
#endif

  libfile = "./" + libfile;

  iface.handle = dlopen(libfile.c_str(), RTLD_NOW);
  if ((iface.handle != NULL) && (dlerror() == NULL))
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
    cerr << "failed" << endl;
#endif

  return success;
}



void
DLLoader::close_library(void* handle)
{
  if (handle != NULL)
    dlclose(handle);
}
