// $Id$

#ifndef _SPECIE_H_
#define _SPECIE_H_

//C++ and boost includes
#include <string>
#include <map>
#include <iostream>
#include <fstream>
#include "boost/assign/list_of.hpp"



//! A class for elemental specie enumeration
/*
 * Methods to overload << >> and = operators are provided
 * to make I/O and comparison operations compatible with
 * strings
 */
class Specie
{
public:

  enum Type {NONE=0, 
                H ,                                                He,
                Li,Be,                              B ,C ,N ,O ,F ,Ne,
                Na,Mg,                              Al,Si,P ,S ,Cl,Ar,
	        K ,Ca,Sc,Ti,V ,Cr,Mn,Fe,Co,Ni,Cu,Zn,Ga,Ge,As,Se,Br,Kr,
                Rb,Sr,Y ,Zr,Nb,Mo,Tc,Ru,Rh,Pd,Ag,Cd,In,Sn,Sb,Te,I ,Xe,
                Cs,Ba,La,Ce,Pr,Nd,Pm,Sm,Eu,Gd,Tb,Dy,Ho,Er,Tm,Yb,
                      Lu,Hf,Ta,W ,Re,Os,Ir,Pt,Au,Hg,Tl,Pb,Bi,Po,At,Rn,
                Fr,Ra,Ac,Th,Pa,U ,Np,Pu,Am,Cm,Bk,Cf,Es,Fm,Md,No,
                      Lr,Rf,Db,Sg,Bh,Hs,Mt,Ds,
                VIRT};


  static
  std::map<Type, std::string> specie_to_string;


   static
   std::map<Type, double> specie_to_mass;


   static
   std::map<std::string, Type> string_to_specie;

    //! Default constructor
    Specie(void);

    //!Default distructor
    ~Specie(void);

    //! Constructor with string specification
    Specie(const std::string& type);

    //! Constructor with enumerator specification
    Specie(const Type& type);

    //!Get mass
    const double get_mass(void) const;

    //!Get specie string
    const std::string& get_string(void) const;

    //! Set specie type as Type enumerator
    void set_type(const Type& type);

    //! Set _type from a string
    void set_type(const std::string& type);

    Specie& operator= (const std::string& type);

    bool operator==(const Specie& specie) const;

    bool operator!=(const Specie& specie) const;


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



inline
const std::string& Specie::get_string(void) const
{
 return Specie::specie_to_string[_type];
}


inline
const double Specie::get_mass(void) const
{
 return Specie::specie_to_mass[_type];
}

inline
void Specie::set_type(const Type& type)
{
 _type = type;
}

inline
void Specie::set_type(const std::string& type)
{
 _type = Specie::string_to_specie[type];
}


inline
bool
Specie::operator==(const Specie& specie) const
{
  return (specie._type == _type);
}

inline
bool
Specie::operator!=(const Specie& specie) const
{
  return !operator==(specie);
}


//Override comparison operator, allows:
// s == 'H'
inline
bool operator== (const Specie& specie, const std::string& type_string)
{
  if (Specie::string_to_specie[type_string] == specie._type) return true;
    else return false;
}

//Override comparison operator, allows:
  // 'H' == s
inline
bool operator== (const std::string& type_string, const Specie& specie)
{
  if (Specie::string_to_specie[type_string] == specie._type) return true;
    else return false;
}

//Override comparison, allow comparison between Specie and Type avoiding
// explicit get_type() call
inline
bool operator== (const Specie& specie, const Specie::Type& type)
{
  return ( specie._type == type);
}

//Override comparison, allow comparison between Specie and Type avoiding
// explicit get_type() call
inline
bool operator== (const Specie::Type& type, const Specie& specie)
{
  return ( specie._type == type);
}

//Override comparison, allow comparison between Specie and Type avoiding
// explicit get_type() call
inline
bool operator!= (const Specie& specie, const Specie::Type& type)
{
  return !(specie == type);
}

//Override comparison, allow comparison between Specie and Type avoiding
// explicit get_type() call
inline
bool operator!= (const Specie::Type& type, const Specie& specie)
{
  return !( type == specie);
}



//Override comparison operator, allows:
  // s == 'H'
inline
bool operator!= (const Specie& specie, const std::string& type_string)
{
  return !(specie == type_string);
}

//Override comparison operator, allows:
  // 'H' == s
inline
bool operator!= (const std::string& type_string, const Specie& specie)
{
  return !(specie == type_string);
}

//Override greater operator, allows:
  // 's1 > s2
inline
bool operator> (const Specie& specie1, const Specie& specie2)
{
  return (specie1._type > specie2._type);
}

//Override greater operator, allows:
  // 's1 < s2
inline
bool operator< (const Specie& specie1, const Specie& specie2)
{
  return (specie1._type < specie2._type);
}


#endif // _SPECIE_H_
