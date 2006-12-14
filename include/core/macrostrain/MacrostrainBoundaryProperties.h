#ifndef _MACROSTRAINBOUNDARYPROPERTIES_H_
#define _MACROSTRAINBOUNDARYPROPERTIES_H_

#include "BoundaryProperties.h"
#include "Material.h"

class MacrostrainBoundaryProperties: public BoundaryProperties
{
 public:

   MacrostrainBoundaryProperties();

 protected:

  virtual void 	do_init (void) = 0;


 private:

  


}; 


//======================================================================//

class MacrostrainSubstrate : public MacrostrainBoundaryProperties
{
 public:
  MacrostrainSubstrate();
  
  ~MacrostrainSubstrate();

  Material* get_material(void) const;

 protected:

  virtual void 	do_init (void);


 private:

  Material*  material;
  

};


inline  Material* MacrostrainSubstrate::get_material( ) const
{

  return material;

}

//=================================================================//

class MacrostrainPressure : public MacrostrainBoundaryProperties
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
