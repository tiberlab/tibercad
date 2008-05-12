// $Id$

#include "TiberEqSystem.h"
#include "TiberNonlinearSystem.h"
#include "InitFailedException.h"

TiberEqSystem::TiberEqSystem(EquationSystems& es,
    const std::string& name, const unsigned int number)
  : ImplicitSystem(es, name, number),
    _type(LINEAR)
{
}



TiberEqSystem*
TiberEqSystem::create(EquationSystems& es,
    const std::string& sysname, const std::string& type,
    const ModelOptions& options)
{
  SystemType sys_type = LINEAR;

  if (type == "linear")
    sys_type = LINEAR;
  else if (type == "nonlinear")
    sys_type = NONLINEAR;
  else
    throw InitFailedException("Unknown equation system type " + sys_type);

  return create(es, sysname, sys_type, options);
}



TiberEqSystem*
TiberEqSystem::create(EquationSystems& es,
    const std::string& sysname, SystemType type,
    const ModelOptions& options)
{
  TiberEqSystem* sys = NULL;

  switch (type)
  {
    case LINEAR:
      break;

    case NONLINEAR:
      sys = TiberNonlinearSystem::create(es, sysname, options);
      break;

    default:
      throw InitFailedException("Unknown equation system type");
  }

  return sys;
}
  
