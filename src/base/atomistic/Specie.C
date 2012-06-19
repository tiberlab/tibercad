#include "Specie.h"

Specie::Specie()
{
  _type = Specie::H;
}

Specie::~Specie(void)
{

}

Specie::Specie(std::string& type)
{
  _type = string_to_specie[type];
}

Specie::Specie(Specie::Type& type)
{
  _type = type;
}

std::map<Specie::Type, std::string>
Specie::specie_to_string = boost::assign::map_list_of
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
       "S");

std::map<std::string, Specie::Type> 
Specie::string_to_specie = boost::assign::map_list_of
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
       ("Ti", Ti) ("W", W) ("U", U) ("V", V) ("Xe", Xe) ("Yb", Yb) ("Y", Y) ("Zn", Zn) ("Zr", Zr) ("S", S);

std::map<Specie::Type, double>
Specie::specie_to_mass = boost::assign::map_list_of
       (Si,28.01) (Ge,72.61) (C,12.0107)(Ga,69.723)(As, 74.92160)(O,15.9994)(Ti, 47.867)(Sr,87.62)(Ca, 40.078)(Al, 26.981538)(N,14.00674)(H,1.00794)(Pb,121.760)(Cl, 35.4527)(Mn,54.938049);

std::map<Specie::Type, double>
Specie::specie_to_mass = boost::assign::map_list_of
       (Si,28.01) (Ge,72.61);


//!Override assignement operator (string input, allows:
//!Specie s = 'H')
Specie& Specie::operator= (const std::string& type)
{
  _type = Specie::string_to_specie[type];
  return *this;
}



//! Override output stream operation (use type as a string)
std::ostream &operator<< (std::ostream& stream,  Specie::Type& type)
{
  stream << Specie::specie_to_string[type];
  return stream;
}

//! Override input stream operation
std::istream &operator>> (std::istream& stream,  Specie::Type& type)
{
  stream >> Specie::specie_to_string[type];
  return stream;
}

//! Override output stream operation (use type as a string)
std::ostream &operator<< (std::ostream& stream, const Specie& specie)
{
  stream << Specie::specie_to_string[specie._type];
  return stream;
}

//! Override input stream operation
std::istream &operator>> (std::istream& stream, Specie& specie)
{
  stream >> Specie::specie_to_string[specie._type];
  return stream;
}
