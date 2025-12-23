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
 * \file TiberModelObject.C
 * \brief tiberCAD API implementation.
 */


#include "tibercad/module/TiberModelObject.h"
#include "tibercad/base/libMeshDefs.h"
#include "tibercad/base/Variable.h"
#include "tibercad/io/Messages.h"
#include "tibercad/base/ModelErrorException.h"

#include "io/DLLoader.h"

#include <vector_value.h>
#include <point.h>

#include <cassert>
#include <sstream>


using namespace std;

USELIBMESHTYPE(RealVectorValue);


TiberModelObject::TiberModelObject(const ModelOptions& options)
  : libMesh::ReferenceCountedObject<TiberModelObject>(),
    _options(options),
    _libhandle(NULL),
    _create(NULL),
    _destroy(NULL),
    _name("")
{
  _name = _options.get_option("name", _name);
}



TiberModelObject::~TiberModelObject(void)
{
  VariableValue::unregister(this);
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
    T& variable, bool override, InitializerBase<T>* initfunc)
{
  string val(_options.get_option(name, ""));
  // if one needs override from strange other sources
  if (override) override_parameter_string(name, val);
  VariableValue::check_and_register(val, variable, this, initfunc);
}


template <typename T>
T
TiberModelObject::get_option(const std::string& name,
    T default_value, bool override) const
{
  string s(_options.get_option(name, ""));
  if (VariableValue::check_string(s))
    throw ModelErrorException("Option \'" + name + "\' cannot "
        " be used as variable (" + s +")");
  if (override) override_parameter_string(name, s);

  T val(default_value);
  if (s.size() > 0) val = Utils::convert<T>(s);

  return val;
}


string
TiberModelObject::get_option(const std::string& name, const char* default_value, bool override) const
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
  if (VariableValue::check_string(s))
    throw ModelErrorException("Option \'" + name + "\' cannot "
        " be used as variable (" + s +")");
  if (override) override_parameter_string(name, s);

  if (s.empty()) return;

  Utils::extract_vector(s, vec);
}


template<>
void
TiberModelObject::get_parameter<RealVectorValue>(const std::string& name,
    RealVectorValue& vec, bool override, InitializerBase<RealVectorValue>*)
{
  string val(_options.get_option(name, ""));
  if (VariableValue::check_string(val))
    throw ModelErrorException("Option \'" + name + "\' cannot "
        " be used as variable (" + val +")");
  // if one needs override from strange other sources
  if (override) override_parameter_string(name, val);

  if (val.empty()) return;

  //VariableValue::check_and_register(val, variable, this, initfunc);
  Utils::extract_vector(val, vec);
}


template<>
void
TiberModelObject::get_parameter<libMesh::RealTensor>(const std::string& name, libMesh::RealTensor& vec, bool override, InitializerBase<libMesh::RealTensor>*)
{
  string val(_options.get_option(name, ""));
  if (VariableValue::check_string(val))
    throw ModelErrorException("Option \'" + name + "\' cannot "
        " be used as variable (" + val +")");
  // if one needs override from strange other sources
  if (override) override_parameter_string(name, val);

  if (val.empty()) return;

  //VariableValue::check_and_register(val, variable, this, initfunc);
  Utils::extract_tensor(val, vec);

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


template <typename T>
void
TiberModelObject::get_option(const std::string& name,
    std::set<T>& vec, bool override) const
{
  vector<T> vvec;
  get_option(name, vvec, override);

  vec.insert(vvec.begin(), vvec.end());
}




void
TiberModelObject::get_option(const std::string& name, libMesh::RealVectorValue& vec, bool override) const
{
  string val(_options.get_option(name, ""));
  if (VariableValue::check_string(val))
    throw ModelErrorException("Option \'" + name + "\' cannot "
        " be used as variable (" + val +")");
  // if one needs override from strange other sources
  if (override) override_parameter_string(name, val);

  if (val.empty()) return;

  //VariableValue::check_and_register(val, variable, this, initfunc);
  Utils::extract_vector(val, vec);
}


void
TiberModelObject::get_option(const std::string& name,
    Point& point, bool override) const
{
  string val(_options.get_option(name, ""));
  if (VariableValue::check_string(val))
    throw ModelErrorException("Option \'" + name + "\' cannot "
        " be used as variable (" + val +")");
  // if one needs override from strange other sources
  if (override) override_parameter_string(name, val);

  if (val.empty()) return;

  //Variable::check_and_register(val, variable, this, initfunc);
  libMesh::RealVectorValue vec(point);
  Utils::extract_vector(val, vec);
  point = vec;
}

void
TiberModelObject::get_option(const std::string& name,
    libMesh::RealTensor& vec, bool override) const
{
  string val(_options.get_option(name, ""));
  if (VariableValue::check_string(val))
    throw ModelErrorException("Option \'" + name + "\' cannot "
        " be used as variable (" + val +")");
  // if one needs override from strange other sources
  if (override) override_parameter_string(name, val);

  if (val.empty()) return;

  //VariableValue::check_and_register(val, variable, this, initfunc);
  Utils::extract_tensor(val, vec);

}



TiberModelObject*
TiberModelObject::_create_from_library(const std::string& name,
    const ModelOptions& options, const void* handle)
{
  TiberModelObject* obj = NULL;

  DLLoader::LibraryInterface iface;

  DLLoader::open_library(name, iface);

  create_t create = (create_t) iface.create_fnc;

  if (create != NULL)
    obj = create(options, handle);

  if (obj != NULL)
  {
    assert(iface.destroy_fnc != NULL);
    assert(iface.handle != NULL);

    obj->_create = create;
    obj->_destroy = (destroy_t) iface.destroy_fnc;
    obj->_libhandle = iface.handle;
  }

  return obj;
}



TiberModelObject*
TiberModelObject::_create_from_object(const TiberModelObject* other, const void* handle)
{
  TiberModelObject* obj = NULL;

  if (other != NULL)
  {
    if (other->_create != NULL)
    obj = other->_create(other->_options, handle);
    if (obj != NULL)
    {
      obj->_create = other->_create;
      obj->_destroy = other->_destroy;
      obj->_libhandle = other->_libhandle;
    }
  }

  return obj;
}



TiberModelObject*
TiberModelObject::create_from_function(create_t create, destroy_t destroy,
    const ModelOptions& options, const void* handle)
{
  assert(create != NULL);
  assert(destroy != NULL);

  TiberModelObject* obj = create(options, handle);

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
  if (p != nullptr)
  {
    Messages m;
    ostringstream os;
    os << "Destroying object: " << p;
    os << ", " << Utils::extract_typename(typeid(*p));
    os << " (" << p->get_name() << ")";
    m.debug(os.str());
    m.indent();

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
    double& val, bool override, InitializerBase<double>* initfunc);

template
void
TiberModelObject::get_parameter<int>(const string& name,
    int& val, bool override, InitializerBase<int>* initfunc);

template
void
TiberModelObject::get_parameter<string>(const string& name,
    string& val, bool override, InitializerBase<string>* initfunc);


template
void
TiberModelObject::get_parameter<unsigned int>(const string& name,
    unsigned int& val, bool override, InitializerBase<unsigned int>* initfunc);

template
void
TiberModelObject::get_parameter<char>(const string& name,
    char& val, bool override, InitializerBase<char>* initfunc);

template
void
TiberModelObject::get_parameter<bool>(const string& name,
    bool& val, bool override, InitializerBase<bool>* initfunc);



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
TiberModelObject::get_option<unsigned short>(const string& name,
    vector<unsigned short>& vec, bool override) const;

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



template
void
TiberModelObject::get_option<double>(const string& name,
    set<double>& vec, bool override) const;

template
void
TiberModelObject::get_option<int>(const string& name,
    set<int>& vec, bool override) const;

template
void
TiberModelObject::get_option<unsigned int>(const string& name,
    set<unsigned int>& vec, bool override) const;

template
void
TiberModelObject::get_option<short>(const string& name,
    set<short>& vec, bool override) const;

template
void
TiberModelObject::get_option<unsigned short>(const string& name,
    set<unsigned short>& vec, bool override) const;


template
void
TiberModelObject::get_option<char>(const string& name,
    set<char>& vec, bool override) const;

template
void
TiberModelObject::get_option<string>(const string& name,
    set<string>& vec, bool override) const;


