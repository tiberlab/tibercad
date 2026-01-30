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
 * \file Specie.C
 * \brief tiberCAD API implementation.
 */

#include "tibercad/atomistic/Specie.h"

Specie::Specie()
{
  _type = NONE;
}

Specie::~Specie(void)
{

}

Specie::Specie(const std::string& type)
{
  (*this) = type;
}

Specie::Specie(const Specie::Type& type)
{
  _type = type;
}

std::map<Specie::Type, std::string>
Specie::_specie_to_string = boost::assign::map_list_of
       (Hf, "Hf") (Al, "Al") (Am, "Am") (Sb, "Sb") (Ag, "Ag") (Ar, "Ar") (As, "As") (Ac, "Ac") (At, "At") (N, "N")
       (Ba, "Ba") (Bk, "Bk") (Be, "Be") (Bi, "Bi") (Bh, "Bh") (B, "B") (Br, "Br") (Cd, "Cd") (Ca, "Ca") (Cf, "Cf")
       (C, "C") (Ce, "Ce") (Cs, "Cs") (Cl, "Cl") (Cr, "Cr") (Co, "Co") (Kr, "Kr") (Cm, "Cm") (Ds, "Ds") (Dy, "Dy")
       (Db, "Db") (Es, "Es") (He, "He") (Er, "Er") (Eu, "Eu") (Fm, "Fm") (Fe, "Fe") (F, "F") (Fr, "Fr") (Gd, "Gd")
       (Ga, "Ga") (Ge, "Ge") (Hs, "Hs") (H, "H") (In, "In") (I, "I") (Ir, "Ir") (La, "La") (Lr, "Lr") (Pb, "Pb")
       (Li, "Li") (Lu, "Lu") (Mg, "Mg") (Mn, "Mn") (Mt, "Mt") (Md, "Md") (Hg, "Hg") (Mo, "Mo") (Nd, "Nd") (Ne,
       "Ne") (Np, "Np") (Ni, "Ni") (Nb, "Nb") (No, "No") (Ho, "Ho") (Au, "Au") (Os, "Os") (O, "O") (Pd, "Pd") (P,
       "P") (Pt, "Pt") (Pu, "Pu") (Po, "Po") (K, "K") (Pr, "Pr") (Pm, "Pm") (Pa, "Pa") (Ra, "Ra") (Rn, "Rn") (Cu,
       "Cu") (Re, "Re") (Rh, "Rh") (Rb, "Rb") (Ru, "Ru") (Rf, "Rf") (Sm, "Sm") (Sc, "Sc") (Sg, "Sg") (Se, "Se")
       (Si, "Si") (Na, "Na") (Sn, "Sn") (Sr, "Sr") (Ta, "Ta") (Tc, "Tc") (Te, "Te") (Tb, "Tb") (Tl, "Tl") (Th, "Th")
       (Tm, "Tm") (Ti, "Ti") (W, "W") (U, "U") (V, "V") (Xe, "Xe") (Yb, "Yb") (Y, "Y") (Zn, "Zn") (Zr, "Zr") (S,
       "S") (VIRT, "Virt") (NONE, "none");

std::map<std::string, Specie::Type> 
Specie::_string_to_specie = boost::assign::map_list_of
       ("Hf", Hf) ("Al", Al) ("Am", Am) ("Sb", Sb) ("Ag", Ag) ("Ar", Ar) ("As", As) ("Ac", Ac) ("At", At) ("N", N)
       ("Ba", Ba) ("Bk", Bk) ("Be", Be) ("Bi", Bi) ("Bh", Bh) ("B", B) ("Br", Br) ("Cd", Cd) ("Ca", Ca) ("Cf", Cf)
       ("C", C) ("Ce", Ce) ("Cs", Cs) ("Cl", Cl) ("Cr", Cr) ("Co", Co) ("Kr", Kr) ("Cm", Cm) ("Ds", Ds) ("Dy", Dy)
       ("Db", Db) ("Es", Es) ("He", He) ("Er", Er) ("Eu", Eu) ("Fm", Fm) ("Fe", Fe) ("F", F) ("Fr", Fr) ("Gd",
       Gd) ("Ga", Ga) ("Ge", Ge) ("Hs", Hs) ("H", H) ("In", In) ("I", I) ("Ir", Ir) ("La", La) ("Lr", Lr) ("Pb", Pb)
       ("Li", Li) ("Lu", Lu) ("Mg", Mg) ("Mn", Mn) ("Mt", Mt) ("Md", Md) ("Hg", Hg) ("Mo", Mo) ("Nd", Nd) ("Ne", Ne)
       ("Np", Np) ("Ni", Ni) ("Nb", Nb) ("No", No) ("Ho", Ho) ("Au", Au) ("Os", Os) ("O", O) ("Pd", Pd) ("P", P)
       ("Pt", Pt) ("Pu", Pu) ("Po", Po) ("K", K) ("Pr", Pr) ("Pm", Pm) ("Pa", Pa) ("Ra", Ra) ("Rn", Rn) ("Cu", Cu)
       ("Re", Re) ("Rh", Rh) ("Rb", Rb) ("Ru", Ru) ("Rf", Rf) ("Sm", Sm) ("Sc", Sc) ("Sg", Sg) ("Se", Se) ("Si", Si)
       ("Na", Na) ("Sn", Sn) ("Sr", Sr) ("Ta", Ta) ("Tc", Tc) ("Te", Te) ("Tb", Tb) ("Tl", Tl) ("Th", Th) ("Tm", Tm)
       ("Ti", Ti) ("W", W) ("U", U) ("V", V) ("Xe", Xe) ("Yb", Yb) ("Y", Y) ("Zn", Zn) ("Zr", Zr) ("S", S)
       ("Virt", VIRT)
       ("none", NONE);

std::map<Specie::Type, double>
Specie::_specie_to_mass = boost::assign::map_list_of
       (Si,28.01) (Ge,72.61) (C,12.0107) (Ga,69.723) (As, 74.92160)
       (O,15.9994)(Ti, 47.867)(Sr,87.62)(Ca, 40.078)(Al, 26.981538)
       (N,14.00674)(H,1.00794)(Pb,207.2)(Cl, 35.4527)(Mn,54.938049)
       (In,114.818)(P,30.973761)(Sb,121.76);


std::vector<double>
Specie::_covalent_radius = { 0.0,
  0.31, 0.28,
  1.28, 0.96, 0.85, 0.76, 0.71, 0.66, 0.57, 0.58,
  1.66, 1.41, 1.21, 1.11, 1.07, 1.05, 1.02, 1.06,
  2.03, 1.76, 1.7, 1.6, 1.53, 1.39, 1.39, 1.32, 1.26, 1.24, 1.32, 1.22, 1.22, 1.2, 1.19, 1.2, 1.2, 1.16,
  2.2, 1.95, 1.9, 1.75, 1.64, 1.54, 1.47, 1.46, 1.42, 1.39, 1.45, 1.44, 1.42, 1.39, 1.39, 1.38, 1.39, 1.4,
  2.44, 2.15, 2.07, 2.04, 2.03, 2.01, 1.99, 1.98, 1.98, 1.96, 1.94, 1.92, 1.92, 1.89, 1.9, 1.87,
  1.87, 1.75, 1.7, 1.62, 1.51, 1.44, 1.41, 1.36, 1.36, 1.32, 1.45, 1.46, 1.48, 1.4, 1.5, 1.5,
  2.6, 2.21, 2.15, 2.06, 2.0, 1.96, 1.9, 1.87, 1.8, 1.69, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0,
  2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0,
  2.0
};

//!Override assignement operator (string input, allows:
//!Specie s = 'H')
Specie& Specie::operator= (const std::string& type)
{
  //_type = _string_to_specie[type];

  if (_string_to_specie.find(type) != _string_to_specie.end())
    _type = _string_to_specie[type];
  else
    _type = VIRT;

  return *this;
}



//! Override output stream operation (use type as a string)
std::ostream &operator<< (std::ostream& stream,  Specie::Type& type)
{
  stream << Specie::_specie_to_string[type];
  return stream;
}

//! Override input stream operation
std::istream &operator>> (std::istream& stream,  Specie::Type& type)
{
  stream >> Specie::_specie_to_string[type];
  return stream;
}

//! Override output stream operation (use type as a string)
std::ostream &operator<< (std::ostream& stream, const Specie& specie)
{
  stream << Specie::_specie_to_string[specie._type];
  return stream;
}

//! Override input stream operation
std::istream &operator>> (std::istream& stream, Specie& specie)
{
  stream >> Specie::_specie_to_string[specie._type];
  return stream;
}
