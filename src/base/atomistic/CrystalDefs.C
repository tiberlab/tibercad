#include "CrystalDefs.h"
#include "Messages.h"
#include "ModelErrorException.h"

#include <fstream>
#include <boost/assign/list_of.hpp>


std::map<std::string, Specie::Type>
CrystalDefs::anion = boost::assign::map_list_of
                                  ("BN", Specie::B)
                                  ("AlN", Specie::N)
                                  ("AlP", Specie::P)
                                  ("AlAs", Specie::As)
                                  ("AlSb", Specie::Sb)
                                  ("GaN", Specie::N)
                                  ("GaP", Specie::P)
                                  ("GaAs", Specie::As)
                                  ("GaSb", Specie::Sb)
                                  ("InN", Specie::N)
                                  ("InP", Specie::P)
                                  ("InAs", Specie::As)
                                  ("InSb", Specie::Sb)
                                  ("SiO2", Specie::O)
                                  ("GeO2", Specie::O)
                                  ("GeS", Specie::S)
                                  ("GeSe", Specie::Se)
                                  ("GeTe", Specie::Te)
                                  ("SnS", Specie::S)
                                  ("SnSe", Specie::Se)
                                  ("SnTe", Specie::Te)
                                  ("PbS", Specie::S)
                                  ("PbSe", Specie::Se)
                                  ("PbTe", Specie::Te)
                                  ("BiTe", Specie::Te)
                                  ("ZnO", Specie::O)
                                  ("ZnS", Specie::S)
                                  ("ZnSe", Specie::Se)
                                  ("ZnTe", Specie::Te)
                                  ("CdO", Specie::O)
                                  ("CdS", Specie::S)
                                  ("CdSe", Specie::Se)
                                  ("CdTe", Specie::Te)
                                  ("HgTe", Specie::Te);


std::map<std::string, Specie::Type>
CrystalDefs::cation = boost::assign::map_list_of
                                  ("AlN", Specie::Al)
                                  ("AlP", Specie::Al)
                                  ("AlAs", Specie::Al)
                                  ("AlSb", Specie::Al)
                                  ("GaN", Specie::Ga)
                                  ("GaP", Specie::Ga)
                                  ("GaAs", Specie::Ga)
                                  ("GaSb", Specie::Ga)
                                  ("InN", Specie::In)
                                  ("InP", Specie::In)
                                  ("InAs", Specie::In)
                                  ("InSb", Specie::In)
                                  ("SiO2", Specie::Si)
                                  ("GeO2", Specie::Ge)
                                  ("GeS", Specie::Ge)
                                  ("GeSe", Specie::Ge)
                                  ("GeTe", Specie::Ge)
                                  ("SnS", Specie::Sn)
                                  ("SnSe", Specie::Sn)
                                  ("SnTe", Specie::Sn)
                                  ("PbS", Specie::Pb)
                                  ("PbSe", Specie::Pb)
                                  ("PbTe", Specie::Pb)
                                  ("BiTe", Specie::Bi)
                                  ("ZnO", Specie::Zn)
                                  ("ZnS", Specie::Zn)
                                  ("ZnSe", Specie::Zn)
                                  ("ZnTe", Specie::Zn)
                                  ("CdO", Specie::Cd)
                                  ("CdS", Specie::Cd)
                                  ("CdSe", Specie::Cd)
                                  ("CdTe", Specie::Cd)
                                  ("HgTe", Specie::Hg);


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


