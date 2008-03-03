#ifndef _MACROSTRAINEXTENDED_H_
#define _MACROSTRAINEXTENDED_H_
//! Substrate boundary condition  
class MacrostrainExtended : public MacrostrainBoundaryProperties
{
 public:
  MacrostrainExtended();
  
 
 

  static MacrostrainExtended* create(void);

 protected:

  virtual void 	do_init (void);

 private:

  
  

};


inline MacrostrainExtended* MacrostrainExtended::create()
{
  return new MacrostrainExtended;
}


#endif
