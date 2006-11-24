#ifndef _ZBSTIFFNESS_H_
#define _ZBSTIFFNESS_H_

#include "Stiffness.h"

class ZbStiffness : public Stiffness
{

 public:

  //!Empty constructor
  ZbStiffness();  

  //!Constructor that sets moduli
  /*! assembles stiffness tensor in crystal system for a zinc-blende crystal 
    \f$ c_{11} = C_{xxxx}, c_{12}=C_{xxyy}, c_{44}=C_{xyxy} \f$ (Voigt notation)
   */ 
  ZbStiffness(double c11, double c12, double c44);

  //!method that sets Young moduli
  void set_moduli(double c11, double c12, double c44);

  //reads database
  void read_database (const Dummy &db);
  
}; 

#endif
