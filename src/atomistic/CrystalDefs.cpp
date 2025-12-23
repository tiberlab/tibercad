/*
 * This file is part of tiberCAD.
 *
 * tiberCAD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * tiberCAD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with tiberCAD. If not, see <https://www.gnu.org/licenses/>.
 */

/*!
 * \file CrystalDefs.C
 * \brief tiberCAD API implementation.
 */

#include "tibercad/atomistic/CrystalDefs.h"
#include "tibercad/io/Messages.h"
#include "tibercad/base/ModelErrorException.h"

#include "boost/assign/list_of.hpp"
#include "boost/lexical_cast.hpp"

#include <fstream>



const std::vector<unsigned int>
CrystalDefs::_class_to_system = 
{
 1, 1,
 2, 2, 2,
 3, 3, 3,
 4, 4, 4, 4, 4, 4, 4,
 5, 5, 5, 5, 5,
 6, 6, 6, 6, 6, 6, 6,
 7, 7, 7, 7, 7,
 0
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
 "23", "m3", "432", "-43m", "m3m",
 ""
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
 "T", "Th", "O", "Td", "Oh",
 ""
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
 
  "tetrahedral", "diploidal", "gyroidal", "hextetrahedral", "hexoctahedral",

  ""
};

const std::vector<std::string>
CrystalDefs::_space_groups =
{
  "P1",

  "P-1",

  "P2", "P2_1", "C2",

  "Pm", "Pc", "Cm", "Cc",

  "P2/m", "P2_1/m", "C2/m", "P2/c", "P2_1/c", "C2/c",

  "P222", "P222_1", "P2_12_12",  "P2_12_12_1", "C222_1",
  "C222", "F222", "I222", "I2_12_12_1",

  "Pmm2", "Pmc2_1", "Pcc2", "Pma2", "Pca2_1", "Pnc2", "Pmn2_1", 
  "Pba2", "Pna2_1", "Pnn2", "Cmm2", "Cmc2_1", "Ccc2", "Amm2",
  "Aem2", "Ama2", "Aea2", "Fmm2", "Fdd2", "Imm2", "Iba2", "Ima2",

  "Pmmm", "Pnnn", "Pccm", "Pban", "Pmma", "Pnna", "Pmna", "Pcca",
  "Pbam", "Pccn", "Pbcm", "Pnnm", "Pmmn", "Pbcn", "Pbca", "Pnma",
  "Cmcm", "Cmce", "Cmmm", "Cccm", "Cmme", "Ccce", "Fmmm", "Fddd",
  "Immm", "Ibam", "Ibca", "Imma",

  "P4", "P4_1", "P4_2", "P4_3", "I4", "I4_1",

  "P-4", "I-4",

  "P4/m", "P4_2/m", "P4/n", "P4_2/n", "I4/m", "I4_1/a",

  "P422", "P42_12", "P4_122", "P4_12_12", "P4_222", "P4_22_12",
  "P4_322", "P4_32_12", "I422", "I4_122",

  "P4mm", "P4bm", "P4_2cm", "P4_2nm", "P4cc", "P4nc", "P4_2mc",
  "P4_2bc", "I4mm", "I4cm", "I4_1md", "I4_1cd",

  "P-42m", "P-42c", "P-42_1m", "P-42_1c", "P-4m2", "P-4c2", "P-4b2",
  "P-4n2", "I-4m2", "I-4c2", "I-42m", "I-42d",

  "P4/mmm", "P4/mcc", "P4/nbm", "P4/nnc", "P4/mbm", "P4/mnc", "P4/nmm",
  "P4/ncc", "P4_2/mmc", "P4_2/mcm", "P4_2/nbc", "P4_2/nnm", "P4_2/mbc",
  "P4_2/mnm", "P4_2/nmc", "P4_2/ncm", "I4/mmm", "I4/mcm", "I4_1/amd", "I4_1/acd",

   "P3", "P3_1", "P3_2", "R3",
   
   "P-3", "R-3",
   
   "P312", "P321", "P3_112", "P3_121", "P3_212", "P3_221", "R32",

   "P3m1", "P31m", "P3c1", "P31c", "R3m", "R3c",
   
   "P-31m", "P-31c", "P-3m1", "P-3c1", "R-3m", "R-3c",
   
   "P6", "P6_1", "P6_5", "P6_2", "P6_4", "P6_3",
   
   "P-6",
   
   "P6/m", "P6_3/m",
   
   "P622", "P6_122", "P6_522", "P6_222", "P6_422", "P6_322",
   
   "P6mm", "P6cc", "P6_3cm", "P6_3mc",
   
   "P-6m2", "P-6c2", "P-62m", "P-62c",
   
   "P6/mmm", "P6/mcc", "P6_3/mcm", "P6_3/mmc",
   
   "P23", "F23", "I23", "P2_13", "I2_13",

   "Pm-3", "Pn-3", "Fm-3", "Fd-3", "Im-3", "Pa-3", "Ia-3",
   
   "P432", "P4_232", "F432", "F4_132", "I432", "P4_332", "P4_132", "I4_132",
   
   "P-43m", "F-43m", "I-43m", "P-43n", "F-43c", "I-43d",
   
   "Pm-3m", "Pn-3n", "Pm-3n", "Pn-3m", "Fm-3m", "Fm-3c",
   "Fd-3m", "Fd-3c", "Im-3m", "Ia-3d",

   "E(3)"
};


std::string
CrystalDefs::convert_to_international_symbol(const std::string& symmetry)
{
  auto is = std::find(_international_sym.begin(),
                      _international_sym.end(), symmetry);
  if (is != _international_sym.end())
    return *is;

  std::string sym = schoenflies_to_IS(symmetry);

  if (sym == "")
    sym = spacegroup_to_IS(symmetry);

  if (sym == "")
    sym = crystal_class_to_IS(symmetry);

  if (sym == "")
  {
    std::string lat = get_bravais_lattice(symmetry);

    if (lat == "aP") sym = "-1";

    if (lat == "mP") sym = "2/m";
    if (lat == "mS") sym = "2/m";

    if (lat == "oP") sym = "mmm";
    if (lat == "oS") sym = "mmm";
    if (lat == "oI") sym = "mmm";
    if (lat == "oF") sym = "mmm";

    if (lat == "tP") sym = "4/mmm";
    if (lat == "tI") sym = "4/mmm";

    if (lat == "hP") sym = "6/mmm";
    if (lat == "hR") sym = "-3m";

    if (lat == "cP") sym = "m3m";
    if (lat == "cI") sym = "m3m";
    if (lat == "cF") sym = "m3m";
  }

  return sym;
}

std::string
CrystalDefs::get_bravais_lattice(const std::string& name)
{
  // her we accept also some sloppy names, for back compatibility
  if ((name == "triclinic") || (name == "aP"))
    return "aP";

  if ((name == "monoclinic") || (name == "monoclinicP") ||
      (name == "mP"))
      return "mP";

  if ((name == "monoclinicS") || (name == "mS"))
      return "mS";

  if ((name == "orthorhombic") || (name == "orthorhombicP") ||
      (name == "oP"))
      return "oP";

  if ((name == "orthorhombicS") || (name == "oS"))
      return "oS";

  if ((name == "orthorhombicI") || (name == "oI"))
      return "oI";

  if ((name == "orthorhombicF") || (name == "oF"))
      return "oF";

  if ((name == "tetragonal") || (name == "tetragonalP") ||
      (name == "tP"))
      return "tP";

  if ((name == "tetragonalI") || (name == "tI"))
      return "tI";

  if ((name == "rhombohedral") || (name == "hR"))
      return "hR";

  if ((name == "hexagonal") || (name == "wz") ||
      (name == "wurtzite") || (name == "hP"))
      return "hP";

  if ((name == "cubic") || (name == "cubicP") ||
      (name == "cP"))
      return "cP";

  if ((name == "bcc") || (name == "cubicI") ||
      (name == "cI"))
      return "cI";

  if ((name == "zb") || (name == "zincblende") ||
      (name == "fcc") || (name == "cubicF") ||
      (name == "cF"))
      return "cF";

  return "";
}



std::string
CrystalDefs::bravais_short_to_long_name(const std::string& in)
{
  std::string name = "";

  switch (in.at(1))
  {
    case 'P':
      name = "primitive ";
      break;

    case 'S':
      name = "base-centered ";
      break;

    case 'I':
      name = "body-centered ";
      break;

    case 'F':
      name = "face-centered ";
      break;

    default:
      break;
  }
  
  switch (in.at(0))
  {
    case 'a':
      name += "triclinic";
      break;

    case 'm':
      name += "monoclinic";
      break;

    case 'o':
      name += "orthorhombic";
      break;

    case 'h':
      if (in.at(1) == 'R')
        name = "rhombohedral";
      else
        name = "hexagonal";
      break;

    case 'c':
      name += "cubic";
      break;

    default:
      break;
  }

  return name;
}



std::string 
CrystalDefs::IS_to_crystal_system(const std::string &int_sym)
{
  std::string system = "";

  auto it = find(_international_sym.begin(), _international_sym.end(), int_sym);
  unsigned int id = std::distance(_international_sym.begin(), it);

  switch (id)
  {
    case 1:
      system = "triclinic";
      break;

    case 2:
      system = "monoclinic";
      break;

    case 3:
      system = "orthorhombic";
      break;

    case 4:
      system = "tetragonal";
      break;

    case 5:
      system = "trigonal";
      break;

    case 6:
      system = "hexagonal";
      break;

    case 7:
      system = "cubic";
      break;

    default:
      break;

  }

  return system;
}



std::string
CrystalDefs::IS_to_crystal_class(const std::string& int_sym)
{
  std::string cc;
  
  auto it = find(_international_sym.begin(), _international_sym.end(), int_sym);
  if (it != _international_sym.end())
  {
    cc = _crystal_classes[std::distance(_international_sym.begin(), it)]; 
  }

  return cc;
}

std::string
CrystalDefs::crystal_class_to_IS(const std::string& cclass)
{
  std::string sym;

  auto it = find(_crystal_classes.begin(), _crystal_classes.end(), cclass);
  if (it != _crystal_classes.end())
  {
    sym = _international_sym[std::distance(_crystal_classes.begin(), it)]; 
  }

  return sym;
}

unsigned int
CrystalDefs::spacegroup_to_IS_id(const std::string &sp_grp)
{
  unsigned int id = 0;
  // sp_grp might be a number
  try
  {
    id = boost::lexical_cast<unsigned int>(sp_grp);
  }
  catch (boost::bad_lexical_cast &)
  {
    auto it = find(_space_groups.begin(), _space_groups.end(), sp_grp);
    id = std::distance(_space_groups.begin(), it);
    
    id++;
  }
  
  unsigned int is_id = 0;

  if (id > 1) ++is_id;
  if (id > 2) ++is_id;
  if (id > 5) ++is_id;
  if (id > 9) ++is_id;
  if (id > 15) ++is_id;
  if (id > 24) ++is_id;
  if (id > 46) ++is_id;
  if (id > 74) ++is_id;
  if (id > 80) ++is_id;
  if (id > 82) ++is_id;
  if (id > 88) ++is_id;
  if (id > 98) ++is_id;
  if (id > 110) ++is_id;
  if (id > 122) ++is_id;
  if (id > 142) ++is_id;
  if (id > 146) ++is_id;
  if (id > 148) ++is_id;
  if (id > 155) ++is_id;
  if (id > 161) ++is_id;
  if (id > 167) ++is_id;
  if (id > 173) ++is_id;
  if (id > 174) ++is_id;
  if (id > 176) ++is_id;
  if (id > 182) ++is_id;
  if (id > 186) ++is_id;
  if (id > 190) ++is_id;
  if (id > 194) ++is_id;
  if (id > 199) ++is_id;
  if (id > 206) ++is_id;
  if (id > 214) ++is_id;
  if (id > 220) ++is_id;
  if (id > 230) ++is_id;

  return is_id;
}


std::string
CrystalDefs::spacegroup_to_IS(const std::string &sp_grp)
{
  unsigned int id =  spacegroup_to_IS_id(sp_grp);

  return _international_sym[id];
}


std::string
CrystalDefs::schoenflies_to_IS(const std::string &schoenflies)
{
  unsigned int id = _international_sym.size() - 1;

  auto it = find(_schoenflies.begin(), _schoenflies.end(), schoenflies);
  if (it != _schoenflies.end())
  {
    id = std::distance(_schoenflies.begin(), it); 
  }
  else
  {
    if (schoenflies == "Ci") id = 1;
    else if (schoenflies == "Cs") id = 3;
    else if (schoenflies == "D2") id = 5;
    else if (schoenflies == "D2h") id = 7;
    else if (schoenflies == "D2d") id = 13;
    else if (schoenflies == "C6i") id = 16;
    else if (schoenflies == "S3") id = 21;
  }

  return _international_sym[id];
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


