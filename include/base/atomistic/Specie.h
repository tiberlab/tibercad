// $Id$

#ifndef _SPECIE_H_
#define _SPECIE_H_

//C++ and boost includes
#include <boost/assign/list_of.hpp>
#include <string>
#include <map>
#include <iostream>
#include <fstream>


enum Type
      {
         Hf, Al, Am, Sb, Ag, Ar, As, Ac, At, N, Ba, Bk, Be, Bi, Bh, B, Br, Cd, Ca, Cf, C,
         Ce, Cs, Cl, Cr, Co, Kr, Cm, Ds, Dy, Db, Es, He, Er, Eu, Fm, Fe, F, Fr, Gd, Ga, Ge,
         Hs, H, In, I, Ir, La, Lr, Pb, Li, Lu, Mg, Mn, Mt, Md, Hg, Mo, Nd, Ne, Np, Ni, Nb,
         No, Ho, Au, Os, O, Pd, P, Pt, Pu, Po, K, Pr, Pm, Pa, Ra, Rn, Cu, Re, Rh, Rb, Ru,
         Rf, Sm, Sc, Sg, Se, Si, Na, Sn, Sr, Ta, Tc, Te, Tb, Tl, Th, Tm, Ti, W, U, V, Xe,
         Yb, Y, Zn, Zr, S
      };

 static
  std::map<Type,std::string> specie_to_string = boost::assign::map_list_of
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


   static
   std::map<std::string, Type> string_to_specie = boost::assign::map_list_of
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

//! A class for elemental specie enumeration
/*
 * Methods to overload << >> and = operators are provided
 * to make I/O and comparison operations compatible with
 * strings
 */

class Specie
{
public:

    //! Default constructor
    Specie(void);

    //!Default distructor
    ~Specie(void);

    //! Constructor with string specification
    Specie(std::string& type);

    //! Constructor with enumerator specification
    Specie(Type& type);

    //! Get reference to specie type
    //const Type& get_type(void) const;

    //!Get specie string
    const std::string& get_string(void) const;

    //! Set specie type as Type enumerator
    void set_type(const Type& type);

    //! Set _type from a string
    void set_type(const std::string& type);

    Specie& operator= (const std::string& type);

    //friend bool operator== (Type& type, std::string& type_string);
    //friend bool operator== (std::string& type_string, Type& type);
    friend bool operator== (const Specie& specie, const std::string& type_string);
    friend bool operator== (const std::string& type_string, const Specie& specie);
    friend bool operator== (const Specie& specie, const Type& type);
    friend bool operator== (const Type& type, const Specie& specie);
    //friend bool operator!= (Type& type, std::string& type_string);
    //friend bool operator!= (std::string& type_string, Type& type);
    friend bool operator!= (const Specie& specie, const Type& type);
    friend bool operator!= (const Type& type, const Specie& specie);
    friend bool operator!= (const Specie& specie, const std::string& type_string);
    friend bool operator!= (const std::string& type_string, const Specie& specie);
    friend bool operator> (const Specie& specie1, const Specie& specie2);
    friend bool operator< (const Specie& specie1, const Specie& specie2);

friend std::ostream& operator<< (std::ostream& stream, Type& type);
friend std::istream& operator>> (std::istream& stream, Type& type);
friend std::ostream& operator<< (std::ostream& stream, const Specie& specie);
friend std::istream& operator>> (std::istream& stream, Specie& specie);
//friend std::ofstream& operator<< (std::ofstream& stream, const Specie& specie);
//friend std::ifstream& operator>> (std::ifstream& stream, Specie& specie);

private:

    Type _type;

};




//inline
//const Type& Specie::get_type(void) const
//{
// return _type;
//}


inline
const std::string& Specie::get_string(void) const
{
 return specie_to_string[_type];
}

inline
void Specie::set_type(const Type& type)
{
 _type = type;
}

inline
void Specie::set_type(const std::string& type)
{
 _type = string_to_specie[type];
}


#endif // _SPECIE_H_
