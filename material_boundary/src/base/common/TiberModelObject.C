// $Id$

#include "TiberModelObject.h"
#include "Variable.h"
#include "DLLoader.h"

#include <cassert>


using namespace std;



TiberModelObject::~TiberModelObject(void)
{
  Variable::unregister(this);
}



bool
TiberModelObject::has_parameter(const std::string& name,
    bool override) const
{
  bool found = _options.find_option(name);

  if (!found && override)
  {
    string s;
    override_parameter_string(name, s);
    if (s.size() > 0)
      found = true;
  }

  return found;
}


template <typename T>
void
TiberModelObject::get_parameter(const std::string& name,
    T& variable, bool override, InitializerBase* initfunc)
{
  string val(_options.get_option(name, ""));
  // if one needs override from strange other sources
  if (override) override_parameter_string(name, val);
  Variable::check_and_register(val, variable, this, initfunc);
}


template <typename T>
T
TiberModelObject::get_option(const std::string& name,
    T default_value, bool override) const
{
  string s(_options.get_option(name, ""));
  if (override) override_parameter_string(name, s);

  T val(default_value);
  if (s.size() > 0) val = Utils::convert<T>(s);

  return val;
}


string
TiberModelObject::get_option(const std::string& name,
    const char* default_value, bool override) const
{
  string val(_options.get_option(name, default_value));
  if (override) override_parameter_string(name, val);

  return val;
}



template <typename T>
void
TiberModelObject::get_parameter(const std::string& name,
    std::vector<T>& vec, bool override)
{
  string s(_options.get_option(name, ""));
  if (override) override_parameter_string(name, s);
  Utils::extract_vector(s, vec);
}



template <typename T>
void
TiberModelObject::get_option(const std::string& name,
    std::vector<T>& vec, bool override) const
{
  string s(_options.get_option(name, ""));
  if (override) override_parameter_string(name, s);
  Utils::extract_vector(s, vec);
}



TiberModelObject*
TiberModelObject::create_from_library(const std::string& name)
{
  TiberModelObject* obj = NULL;

  DLLoader::LibraryInterface iface;

  bool success = DLLoader::open_library(name, iface);

  if (success)
  {
    create_t create = (create_t) iface.create_fnc;

    assert(create != NULL);
    obj = create();

    if (obj != NULL)
    {
      assert(iface.destroy_fnc != NULL);
      assert(iface.handle != NULL);

      obj->_create = create;
      obj->_destroy = (destroy_t) iface.destroy_fnc;
      obj->_libhandle = iface.handle;
    }
  }

  return obj;
}


TiberModelObject*
TiberModelObject::create_from_function(create_t create, destroy_t destroy)
{
  assert(create != NULL);
  assert(destroy != NULL);

  TiberModelObject* obj = create();

  if (obj != NULL)
  {
    obj->_create = create;
    obj->_destroy = destroy;
  }

  return obj;
}


void
TiberModelObject::destroy(TiberModelObject* p)
{
  if (p != NULL)
  {

    libhandle_t libhandle = p->_libhandle;
    destroy_t destroy_fnc = p->_destroy;

    if (destroy_fnc != NULL)
      destroy_fnc(p);
    else
      delete p;

    if (libhandle != NULL)
      DLLoader::close_library(libhandle);
  }
}


//
// explicit instantiations
//

template
void
TiberModelObject::get_parameter<double>(const string& name,
    double& val, bool override, InitializerBase* initfunc);

template
void
TiberModelObject::get_parameter<int>(const string& name,
    int& val, bool override, InitializerBase* initfunc);

template
void
TiberModelObject::get_parameter<string>(const string& name,
    string& val, bool override, InitializerBase* initfunc);


template
void
TiberModelObject::get_parameter<unsigned int>(const string& name,
    unsigned int& val, bool override, InitializerBase* initfunc);

template
void
TiberModelObject::get_parameter<char>(const string& name,
    char& val, bool override, InitializerBase* initfunc);

template
void
TiberModelObject::get_parameter<bool>(const string& name,
    bool& val, bool override, InitializerBase* initfunc);



template
double
TiberModelObject::get_option<double>(const string& name,
    double val, bool override) const;

template
int
TiberModelObject::get_option<int>(const string& name,
    int val, bool override) const;

template
unsigned int
TiberModelObject::get_option<unsigned int>(const string& name,
    unsigned int val, bool override) const;

template
short
TiberModelObject::get_option<short>(const string& name,
    short val, bool override) const;


template
bool
TiberModelObject::get_option<bool>(const string& name,
    bool val, bool override) const;

template
char
TiberModelObject::get_option<char>(const string& name,
    char val, bool override) const;

template
string
TiberModelObject::get_option<string>(const string& name,
    string val, bool override) const;

template
const char*
TiberModelObject::get_option<const char*>(const string& name,
    const char* val, bool override) const;




template
void
TiberModelObject::get_option<double>(const string& name,
    vector<double>& vec, bool override) const;

template
void
TiberModelObject::get_option<int>(const string& name,
    vector<int>& vec, bool override) const;

template
void
TiberModelObject::get_option<unsigned int>(const string& name,
    vector<unsigned int>& vec, bool override) const;

template
void
TiberModelObject::get_option<short>(const string& name,
    vector<short>& vec, bool override) const;

template
void
TiberModelObject::get_option<bool>(const string& name,
    vector<bool>& vec, bool override) const;

template
void
TiberModelObject::get_option<char>(const string& name,
    vector<char>& vec, bool override) const;

template
void
TiberModelObject::get_option<string>(const string& name,
    vector<string>& vec, bool override) const;




template
void
TiberModelObject::get_parameter<double>(const string& name,
    vector<double>& vec, bool override);

template
void
TiberModelObject::get_parameter<int>(const string& name,
    vector<int>& vec, bool override);

template
void
TiberModelObject::get_parameter<unsigned int>(const string& name,
    vector<unsigned int>& vec, bool override);

template
void
TiberModelObject::get_parameter<short>(const string& name,
    vector<short>& vec, bool override);

template
void
TiberModelObject::get_parameter<bool>(const string& name,
    vector<bool>& vec, bool override);

template
void
TiberModelObject::get_parameter<char>(const string& name,
    vector<char>& vec, bool override);

template
void
TiberModelObject::get_parameter<string>(const string& name,
    vector<string>& vec, bool override);
