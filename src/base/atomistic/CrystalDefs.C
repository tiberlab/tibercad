#include "CrystalDefs.h"
#include "Messages.h"
#include "ModelErrorException.h"

#include <fstream>
#include <boost/assign/list_of.hpp>



const std::vector<unsigned int>
CrystalDefs::_class_to_system = 
{
 1, 1,
 2, 2, 2,
 3, 3, 3,
 4, 4, 4, 4, 4, 4, 4,
 5, 5, 5, 5, 5,
 6, 6, 6, 6, 6, 6, 6,
 7, 7, 7, 7, 7
};

const std::vector<std::string>
CrystalDefs::_international_sym =
{
 "1", "-1",
 "2", "m", "2/m",
 "222", "mm2", "mmm",
 "4", "-4", "4/m", "422", "4mm", "-42m", "4/mmm",
 "3", "-3", "32", "3m", "-3m",
 "6", "-6", "6/m", "622", "6mm", "-6m2", "6/mmm",
 "23", "m3", "432", "-43m", "m3m"
};


const std::vector<std::string>
CrystalDefs::_schoenflies =
{
 "C1", "S2",
 "C2", "C1h", "C2h",
 "V", "C2v", "Vh",
 "C4", "S4", "C4h", "D4", "C4v", "Vd", "D4h",
 "C3", "S6", "D3", "C3v", "D3d",
 "C6", "C3h", "C6h", "D6", "C6v", "D3h", "D6h",
 "T", "Th", "O", "Td", "Oh"
};

const std::vector<std::string>
CrystalDefs::_crystal_classes =
{
  "triclinic-pedial", "triclinic-pinacoidal",

  "monoclinic-sphenoidal", "monoclinic-domatic", "monoclinic-prismatic",

  "orthorhombic-disphenoidal", "orthorhombic-pyramidal", "orthorhombic-dipyramidal",

  "tetragonal-pyramidal", "tetragonal-disphenoidal", "tetragonal-dipyramidal",
  "tetragonal-trapezoidal", "ditetragonal-pyramidal", "tetragonal-scalenoidal",
  "ditetragonal-dipyramidal",

  "trigonal-pyramidal", "rhombohedral", "trigonal-trapezoidal", "ditrigonal-pyramidal",
  "ditrigonal-scalahedral",

  "hexagonal-pyramidal", "trigonal-dipyramidal", "hexagonal-dipyramidal",
  "hexagonal-trapezoidal", "dihexagonal-pyramidal", "ditrigonal-dipyramidal",
  "dihexagonal-dipyramidal",
 
  "tetrahedral", "diploidal", "gyroidal", "hextetrahedral", "hexoctahedral"
};


std::string
CrystalDefs::IS_to_crystal_class(const std::string& int_sym)
{
  std::string cc;
  
  auto it = find(_international_sym.begin(), _international_sym.end(), int_sym);
  if (it != _international_sym.end())
  {
    cc = _crystal_classes[std::distance(it, _international_sym.begin())]; 
  }

  return cc;
}


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


