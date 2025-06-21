// $Id: TmmBulkModel.C 4391 2017-04-07 11:16:58Z pecchia $

#include "TmmDipoleSource.h"
#include "Database.h"
#include "Messages.h"

#include <boost/filesystem/operations.hpp>

using std::string;
using namespace libMesh;


TmmDipoleSource::TmmDipoleSource(const ModelOptions& options) :
  PhysicalModel(options),
  _emission_power(0.0)
{


}


const double&
TmmDipoleSource::get_emission_power(void) const
{

  return(_emission_power);
}


void
TmmDipoleSource::set_emission_power(const double& emission_power)
{
  _emission_power = emission_power;
}




