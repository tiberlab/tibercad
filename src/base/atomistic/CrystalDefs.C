#include "CrystalDefs.h"
#include "Messages.h"
#include "ModelErrorException.h"

#include <fstream>
#include <boost/assign/list_of.hpp>


std::map<std::string, Specie::Type>
CrystalDefs::anion = boost::assign::map_list_of("AlAs", Specie::As)
                                  ("AlN", Specie::N)
                                  ("AlP", Specie::P)
                                  ("AlSb", Specie::Sb)
                                  ("GaAs", Specie::As)
                                  ("GaN", Specie::N)
                                  ("GaP", Specie::P)
                                  ("GaSb", Specie::Sb)
                                  ("InAs", Specie::As)
                                  ("InN", Specie::N)
                                  ("InP", Specie::P)
                                  ("InSb", Specie::Sb);

std::map<std::string, Specie::Type>
CrystalDefs::cation = boost::assign::map_list_of("AlAs", Specie::Al)
                                  ("AlN", Specie::Al)
                                  ("AlP", Specie::Al)
                                  ("AlSb", Specie::Al)
                                  ("GaAs", Specie::Ga)
                                  ("GaN", Specie::Ga)
                                  ("GaP", Specie::Ga)
                                  ("GaSb", Specie::Ga)
                                  ("InAs", Specie::In)
                                  ("InN", Specie::In)
                                  ("InP", Specie::In)
                                  ("InSb", Specie::In);


bool CrystalDefs::is_anion(const std::string mat_name, const Specie sp)
{
  std::map<const std::string, Specie::Type>::iterator it(
      CrystalDefs::anion.find(mat_name));

  if (it == CrystalDefs::anion.end())
  {
    std::ostringstream os;
    os << "CrystalDefs: Anion not defined for material ";
    os << mat_name;
    throw(ModelErrorException(os.str()));
  }

  return(it->second == sp);
}


bool CrystalDefs::is_cation(const std::string mat_name, const Specie sp)
{
  std::map<const std::string, Specie::Type>::iterator it(
      CrystalDefs::cation.find(mat_name));

  if (it == CrystalDefs::cation.end())
  {
    std::ostringstream os;
    os << "CrystalDefs: Cation not defined for material ";
    os << mat_name;
    throw(ModelErrorException(os.str()));
  }

  return(it->second == sp);
}


