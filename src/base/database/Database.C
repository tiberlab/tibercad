// $Id$

#include "Database.h"


bool
Database::is_alloy(const std::string& name) const
{
  if (name == "AlGaAs")  return true;
  if (name == "InGaAs")  return true;
  if (name == "AlGaN")  return true;
  if (name == "InGaN")  return true;
  if (name == "AlInN")  return true;

  return false;
}

void
Database::get_alloy_components(const std::string& alloy,
    std::string& comp_A, std::string& comp_B) const
{
  if (alloy == "AlGaAs") { comp_A = "AlAs"; comp_B = "GaAs"; }
  if (alloy == "InGaAs") { comp_A = "InAs"; comp_B = "GaAs"; }
  if (alloy == "AlGaN") { comp_A = "AlN"; comp_B = "GaN"; }
  if (alloy == "InGaN") { comp_A = "InN"; comp_B = "GaN"; }
  if (alloy == "AlInN") { comp_A = "AlN"; comp_B = "InN"; }
}
