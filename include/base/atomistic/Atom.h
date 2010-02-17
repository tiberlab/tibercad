#ifndef _ATOM_H_
#define _ATOM_H_

#include "tensor.h"
#include "TypeDefs.h"
#include "point.h"
#include "dof_object.h"

//C++ and boost includes
#include <boost/assign/list_of.hpp>

class Elem;


//Enumerator for atomic specie and maps for getting string from enumeraqtor and enumerator from string
enum Specie
{
   Hf, Al, Am, Sb, Ag, Ar, As, Ac, At, N, Ba, Bk, Be, Bi, Bh, B, Br, Cd, Ca, Cf, C,
   Ce, Cs, Cl, Cr, Co, Kr, Cm, Ds, Dy, Db, Es, He, Er, Eu, Fm, Fe, F, Fr, Gd, Ga, Ge,
   Hs, H, In, I, Ir, La, Lr, Pb, Li, Lu, Mg, Mn, Mt, Md, Hg, Mo, Nd, Ne, Np, Ni, Nb,
   No, Ho, Au, Os, O, Pd, P, Pt, Pu, Po, K, Pr, Pm, Pa, Ra, Rn, Cu, Re, Rh, Rb, Ru,
   Rf, Sm, Sc, Sg, Se, Si, Na, Sn, Sr, Ta, Tc, Te, Tb, Tl, Th, Tm, Ti, W, U, V, Xe,
   Yb, Y, Zn, Zr, S
};


static
std::map<Specie,std::string> specie_to_string = boost::assign::map_list_of
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
std::map<std::string, Specie> string_to_specie = boost::assign::map_list_of
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


//! Contains Atom definition
/*!
 * Atom is defined mainly by atomic specie and
 * spatial vector giving the position (from library
 * tensor.h)
 */
class Atom
{
public:

  //! Atom constructor
  Atom();

  //! Constructor with specie and position initializations
  Atom(std::string& init_specie, Tensor1& init_position);

  //!Complete constructor: specifies all atom characteristics
  Atom(std::string& specie, Tensor1& position, unsigned int flag);

  //! Atom destructor
  ~Atom();

  //! Set atom specie name
  void set_specie(const std::string& sp);

  //! Get atom specie name
  const std::string& get_specie() const;

  //!Set atom position (1X3 Tensor is used)
  void set_position(const Tensor1 pos);

  //! Get position coordinate i (x=1, y=2, z=3)
  double get_position(int i) const;

  //! Get the whole position (1X3 Tensor)
  Tensor1 get_position() const;

  //! Get the ID of the region the atom belongs to
  //! (Note: little errors may occur using basis or conventional cell
  //! preservation in Atomistic Generator)
  int get_region_ID() const;

  //! Set a general purpose integer flag, used internally. 0 is default safe value
  void set_flag(const unsigned int fg);

  //! Get the general purpose flag
  unsigned int get_flag() const;

  //! True if atom belong to structure. Useful during structure construction
  bool belong_to_structure;

  //!Set element
  void set_elem(Elem* el);

  //!Get element
  const Elem* get_elem() const;

private:

  //!Element atom belongs to
  Elem* _el;

  //! Atomic specie (short name)
  std::string _specie;

  //! Atom position
  Point _position;

  //! A general purpose integer flag (for example used in passivation)
  unsigned int _flag;

};


//----------------------------------------------------
// Inline member functions
//----------------------------------------------------

inline
void Atom::set_specie(const std::string& sp)
{
  _specie=sp;
}


inline
const std::string& Atom::get_specie(void) const
{
  return _specie;
}


inline
double Atom::get_position(int i) const
{
  return _position(i - 1);
}


inline
void Atom::set_flag(const unsigned int fg)
{
  _flag = fg;
}


inline
unsigned int Atom::get_flag() const
{
  return _flag;
}


inline
const Elem* Atom::get_elem() const
{
  return _el;
}

#endif // _ATOM_H_

