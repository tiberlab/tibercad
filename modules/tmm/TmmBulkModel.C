/*  
 * This file is part of the tiberCAD module tmm.
 *
 * tiberCAD modules are licensed under the GNU General Public License v3.
 *
 * tiberCAD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published
 * by the Free Software Foundation, either version 3 of the License,
 * or (at your option) any later version.
 *
 * tiberCAD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with tiberCAD. If not, see <http://www.gnu.org/licenses/>.
 */

/*!
 * \file TmmBulkModel.C
 * \brief tiberCAD tmm module implementation.
 *
 * \note This file is part of module tmm.
 */


#include "TmmBulkModel.h"
#include "tibercad/io/Database.h"
#include "tibercad/io/Messages.h"

#include <boost/filesystem/operations.hpp>

using std::string;
using namespace libMesh;


TiberModelObject*
TmmBulkModel::_create(const ModelOptions& options, const void*)
{
  return new TmmBulkModel(options);
}


inline
TmmBulkModel::TmmBulkModel(const ModelOptions& options) :
  PhysicalModel(options)
{


}



void
TmmBulkModel::_destroy(TiberModelObject* p)
{
  delete p;
}


TmmBulkModel*
TmmBulkModel::create(const Material* mat, const ModelOptions& options)
{
  string type("default");
  options.get_option("type", type);

  TmmBulkModel* pm = NULL;
  if (type == "default")
  {
    // we create the default model from explicit creation method
    pm = PhysicalModel::create<TmmBulkModel>(_create, _destroy, mat, options);

  }
  else
  {
    // there is no such model, at the moment
    type = "bulk_model" + type;
    pm = PhysicalModel::create<TmmBulkModel>(type, mat, options);
  }

  return(pm);

}

std::pair<double, double>
TmmBulkModel::interpolate(double wavelength) const
{
  // calculate index for adressing array
  double dx = (_wavelengths.back() - _wavelengths.front()) / _addressing.size();
  int index = std::max(0, static_cast<int>(floor((wavelength - _wavelengths.front()) / dx)));
  index = std::min(index, static_cast<int>(_addressing.size() - 1));

  unsigned int ctr = _addressing[index];

  while ((ctr < _wavelengths.size()) && (_wavelengths[ctr] < wavelength))
    ++ctr;

  double n, k;

  if (ctr == 0)
  {
    n = _n_data[0];
    k = _k_data[0];
  }
  else if (ctr == _wavelengths.size())
  {
    n = _n_data.back();
    k = _k_data.back();
  }
  else
  {
    double x1 = _wavelengths[ctr - 1];
    double x2 = _wavelengths[ctr];

    double frac = (wavelength - x1) / (x2 - x1);
    n = (_n_data[ctr] - _n_data[ctr - 1]) * frac + _n_data[ctr - 1];
    k = (_k_data[ctr] - _k_data[ctr - 1]) * frac + _k_data[ctr - 1];
  }

  return(std::make_pair(n, k));
}


libMesh::Complex
TmmBulkModel::get_permittivity(double lambda) const
{
  auto nk = interpolate(lambda);

  libMesh::Complex eps(nk.first, nk.second);
  eps = eps * eps;

  return(eps);
}

libMesh::Complex
TmmBulkModel::get_refractive_index(double lambda) const
{
  auto nk = interpolate(lambda);

  libMesh::Complex ri(nk.first, nk.second);

  return(ri);
}

double
TmmBulkModel::get_coherent_index(void) const
{

  double ri = _incoherent_index;

  return(ri);
}

double
TmmBulkModel::get_emission_power(void) const
{

  double ri = _emission_power;

  return(ri);
}



void
TmmBulkModel::read_database(void)
{
  const Database& db = get_database();
  db.set_section("permittivity");

  string source = db.get("optical_data", "");

  // open the given file
  // for this we assume that optical_data is a relative path, relative to the database file
  string dbfile = db.get_data_file();
  boost::filesystem::path p(dbfile);
  p = p.parent_path() / boost::filesystem::path(source);

  if (!boost::filesystem::exists(p) || !boost::filesystem::is_regular_file(p))
  {
    std::ostringstream os;
    os << "Material database for " << db.get_material() << " does not provide optical data";
    throw InitFailedException(os.str());
  }

  _datafile = p.string();


  // now read the file
  // structure is:
  // # comments
  // # comments
  // wavelength n k

  // a line buffer
  const size_t buf_len = 256;
  char buf[buf_len];

  std::ifstream is(_datafile);
  if (is.fail() || !is.good())
    throw InitFailedException("Cannot read data from \"" + _datafile + "\"");

  size_t i = 0;
  while (is.good())
  {
    if (i == _wavelengths.size())
    {
      size_t n_new = _wavelengths.size() + 100;
      _wavelengths.reserve(n_new);
      _n_data.reserve(n_new);
      _k_data.reserve(n_new);
    }

    is.getline(buf, buf_len);
    if ((buf[0] != '#') && (buf[0] != '%') && (buf[0] != '/'))
    {
      std::istringstream in(buf);

      double w, n, k;
      if (in >> w)
      {
        if ((in >> n) && (in >> k))
        {
          _wavelengths.push_back(w);
          _n_data.push_back(n);
          _k_data.push_back(k);
          i++;
        }
      }
    }
  }
  is.close();

  _wavelengths.resize(_wavelengths.size());
  _n_data.resize(_n_data.size());
  _k_data.resize(_k_data.size());

  // _addressing is used for faster access during interpolation
  double min = _wavelengths.front();
  double max = _wavelengths.back();
  const int N = 50;
  double dx = (max - min) / N;
  _addressing.resize(N);

  unsigned int ctr = 0;
  for (unsigned int i = 0; i < N; ++i)
  {
    double x = min + i * dx;

    while ((_wavelengths[ctr] < x) && (ctr < _wavelengths.size()))
      ++ctr;

    _addressing[i] = ctr;
  }
}


void
TmmBulkModel::do_init(void)
{

  //SubmodelIterator it = submodels_begin("InCoherent");
  //InCoherentModel* pm =  dynamic_cast<InCoherentModel*>(it->second);
  //_incoherent_index = pm->get_InCoherent_Index();
  get_parameter("Incoherency", _incoherent_index);
  // string gen_str(get_option("generation", "0"));
  // istringstream is(gen_str);
  // std::cout<<"gen_str is :"<<gen_str<<std::endl;
  // std::cout<<"is is :"<<is<<std::endl;

}


void
TmmBulkModel::do_print_info(void)
{
  Messages::info("Reading optical data from " + _datafile);
  std::ostringstream os;
  os << "wavelength range: " << _wavelengths.front() << " - " << _wavelengths.back() << " nm";
  Messages::info(os.str());
}


void
TmmBulkModel::prepare_submodels(void)
{

  create_submodels(_DS, "dipole_source");
}


void
TmmBulkModel::calculate(const libMesh::Elem* elem, const libMesh::Point& point, double lambda)
{
  _emission_power = 0;
  for (ID n = 0; n < _DS.size(); n++)
  {
    _DS[n]->calculate(elem, point, lambda);
    _emission_power += _DS[n]->get_emission_power();
  }

}
