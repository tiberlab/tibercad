// $Id$

#include "tibercad/profiles/ExternalProfile.h"
#include "tibercad/profiles/ExtProfile1D.h"
#include "tibercad/profiles/UniformRandomAlloy.h"
#include "tibercad/profiles/LinearProfile.h"
#include "tibercad/profiles/ExponentialProfile.h"
#include "tibercad/profiles/GaussianProfile.h"
#include "tibercad/profiles/CompositeProfile.h"
#include "tibercad/module/SimulationInterface.h"

using namespace std;


ExternalProfile::ExternalProfile(const ModelOptions& options) :
  TiberModelObject(options)
{

}


ExternalProfile::~ExternalProfile(void)
{
}

void
ExternalProfile::setup(void)
{
  string src = get_option("data_source", "");
  _data_source = SimulationInterface::find_solution_provider(src);
  if (_data_source.second == INVALID_ID)
    throw InitFailedException("'" + src + "' is invalid data source for external profile");
}

ExternalProfile*
ExternalProfile::create(const ModelOptions& options)
{
  ExternalProfile* pr = nullptr;

  if (options.get_name() == "1d_profile")
    pr = new ExtProfile1D(options);
  else if (options.get_name() == "uniform_random")
    pr = new UniformRandomAlloy(options);
  else if (options.get_name() == "linear")
    pr = new LinearProfile(options);
  else if (options.get_name() == "exponential")
    pr = new ExponentialProfile(options);
  else if (options.get_name() == "gaussian")
    pr = new GaussianProfile(options);
  else if (options.get_name() == "composite")
    pr = new CompositeProfile(options);
  else
  {
    pr = new ExternalProfile(options);
    pr->setup();
  }

  return pr;
}


double
ExternalProfile::get_data(const Elem* elem) const
{
  double conc = 0.0;
  int n_nodes = elem->n_nodes();

  for (int i = 0; i < n_nodes; ++i)
    conc += get_data(elem, elem->point(i));


  return(conc / n_nodes);
}

double
ExternalProfile::get_data(const Elem* elem, const Point& p) const
{
  // here we arrive only if data source is valid
  vector<double> tmp(1);
  _data_source.first->get_solution(elem, _data_source.second, tmp, vector<Point>(1, p));

  return(tmp[0]);
}

std::pair<double, double>
ExternalProfile::get_min_max(void) const
{
  return(make_pair(numeric_limits<double>::min(), numeric_limits<double>::max()));
}
