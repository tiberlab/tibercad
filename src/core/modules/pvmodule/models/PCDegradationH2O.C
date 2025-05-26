#include "PCDegradationH2O.h"

#include "TiberModule.h"

PCDegradationH2O::PCDegradationH2O(const ModelOptions& options)
  : Photocurrent(options)
{
}

PCDegradationH2O*
PCDegradationH2O::create(const ModelOptions& options)
{
  return new PCDegradationH2O(options);
}

void
PCDegradationH2O::do_init(void)
{
  _initial_current = get_option("initial_photocurrent", _initial_current);
}


double
PCDegradationH2O::do_get_photocurrent(const libMesh::Elem* elem,
                                      const libMesh::Point& p) const
{
  return _initial_current;
}
