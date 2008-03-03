#ifndef _MACROSTRAINSUBSTRATE_H_
#define _MACROSTRAINSUBSTRATE_H_
//! Substrate boundary condition  
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
//! returns a pointer to the substrate material 
inline  Material* MacrostrainSubstrate::get_material( ) const
{

  return material;

}

#endif
