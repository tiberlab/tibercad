#ifndef _MACROSTRAINBOUNDARYPROPERTIES_H_
#define _MACROSTRAINBOUNDARYPROPERTIES_H_

#include "BoundaryProperties.h"
#include "Material.h"

class MacrostrainBoundaryProperties: public BoundaryProperties
{
 public:

   MacrostrainBoundaryProperties();

   static  MacrostrainBoundaryProperties* create(const std::string & name,  const ModelOptions &   options );

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

  static MacrostrainSubstrate* create(void);

 protected:

  virtual void 	do_init (void);


 private:

  Material*  material;
  

};


inline MacrostrainSubstrate* MacrostrainSubstrate::create()
{
  return new MacrostrainSubstrate;
}

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

  static MacrostrainPressure* create();

 protected:

  virtual void 	do_init (void);

  


 private:

  //! value of pressure in GPa
  double value;

}; 

inline MacrostrainPressure* MacrostrainPressure::create()
{
  return new MacrostrainPressure;
}


inline double MacrostrainPressure::get_value(void) const
{
  return value;
}


//================================================================//
#endif
