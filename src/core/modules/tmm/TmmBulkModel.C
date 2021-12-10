// $Id: TmmBulkModel.C 4391 2017-04-07 11:16:58Z pecchia $

#include "TmmBulkModel.h"
#include "Database.h"
#include "Messages.h"
//#include "PermittivityModel.h"

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

  return(std::make_pair(n,k));
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
TmmBulkModel::get_incoherance_index() const
{
  return(0);
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
  // Maybe it would be more elegant to extend the existing
  // permittivity model implementation
  //ModelOptions opts;
  //opts.set_option("type", "constant");
  //create_submodel(_permittivity_model, "permittivity", opts);
  
  // alternative way to create internal submodels:
  //
  // PermittivityModel* mod = PhysicalModel::create("permittivity", opts);
  // add_submodel("permittivity", mod)

  // NOTE: all submodels are initialized automatically before calling do_init()
}
