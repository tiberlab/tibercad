#include "Specie.h"

Specie::Specie()
{
  _type = H;
}

Specie::~Specie(void)
{

}

Specie::Specie(std::string& type)
{
  _type = string_to_specie[type];
}

Specie::Specie(Type& type)
{
  _type = type;
}


//!Override assignement operator (string input, allows:
//!Specie s = 'H')
Specie& Specie::operator= (const std::string& type)
{
  _type = string_to_specie[type];
  return *this;
}


//!Override comparison operator, allows:
//! s == 'H'
bool operator== (const Specie& specie, const std::string& type_string)
{
  if (string_to_specie[type_string] == specie._type) return true;
    else return false;
}

//!Override comparison operator, allows:
  //! 'H' == s
bool operator== (const std::string& type_string, const Specie& specie)
{
  if (string_to_specie[type_string] == specie._type) return true;
    else return false;
}

//!Override comparison, allow comparison between Specie and Type avoiding
//! explicit get_type() call
bool operator== (const Specie& specie, const Type& type)
    {
  return ( specie._type == type);
    }

//!Override comparison, allow comparison between Specie and Type avoiding
//! explicit get_type() call
bool operator== (const Type& type, const Specie& specie)
    {
  return ( specie._type == type);
    }

//!Override comparison, allow comparison between Specie and Type avoiding
//! explicit get_type() call
bool operator!= (const Specie& specie, const Type& type)
    {
  return !(specie == type);
    }

//!Override comparison, allow comparison between Specie and Type avoiding
//! explicit get_type() call
bool operator!= (const Type& type, const Specie& specie)
    {
  return !( type == specie);
    }

////!Override comparison operator, allows:
//    //! s.get_type() == 'H'
//bool operator!= (Type& type, std::string& type_string)
//{
//  return !(type == string_to_specie[type_string]);
//}
//
////!Override comparison operator, allows:
//    //! 'H' == s.get_type()
//bool operator!= (std::string& type_string, Type& type)
//{
//  return !(type == type_string);
//}

//!Override comparison operator, allows:
  //! s == 'H'
bool operator!= (const Specie& specie, const std::string& type_string)
{
  return !(specie == type_string);
}

//!Override comparison operator, allows:
  //! 'H' == s
bool operator!= (const std::string& type_string, const Specie& specie)
{
  return !(specie == type_string);
}

//!Override greater operator, allows:
  //! 's1 > s2
bool operator> (const Specie& specie1, const Specie& specie2)
{
  return (specie1._type > specie2._type);
}

//!Override greater operator, allows:
  //! 's1 < s2
bool operator< (const Specie& specie1, const Specie& specie2)
{
  return (specie1._type < specie2._type);
}

//! Override output stream operation (use type as a string)
std::ostream &operator<< (std::ostream& stream, Type& type)
{
  stream << specie_to_string[type];
  return stream;
}

//! Override input stream operation
std::istream &operator>> (std::istream& stream, Type& type)
{
  stream >> specie_to_string[type];
  return stream;
}

//! Override output stream operation (use type as a string)
std::ostream &operator<< (std::ostream& stream, const Specie& specie)
{
  stream << specie_to_string[specie._type];
  return stream;
}

//! Override input stream operation
std::istream &operator>> (std::istream& stream, Specie& specie)
{
  stream >> specie_to_string[specie._type];
  return stream;
}
