#ifndef _MACROSTRAINBOUNDARYPROPERTIES_H_
#define _MACROSTRAINBOUNDARYPROPERTIES_H_

#include "BoundaryProperties.h"
#include "RotatedCrystal.h"

class MacrostrainBoundaryProperties: public BoundaryProperties
{
 public:

   MacrostrainBoundaryProperties();

 protected:

  virtual void 	do_init (void) = 0;


 private:

  


}; 


//======================================================================//

class MacrostrainSubstrate : MacrostrainBoundaryProperties
{
 public:
  MacrostrainSubstrate();
  
  ~MacrostrainSubstrate();

  RotatedCrystal* get_crystal(void) const;

 protected:

  virtual void 	do_init (void);


 private:
  RotatedCrystal* crystal;
  

};


inline  RotatedCrystal* MacrostrainSubstrate::get_crystal( ) const
{

  return crystal;

}

//=================================================================//

class MacrostrainPressure : MacrostrainBoundaryProperties
{
 public:

  
  MacrostrainPressure();

  double get_value(void) const;

 protected:

  virtual void 	do_init (void);

  


 private:

  //! value of pressure in GPa
  double value;

}; 



inline double MacrostrainPressure::get_value(void) const
{
  return value;
}


//================================================================//
#endif
