#ifndef _MACROSTRAINEXTENDED_H_
#define _MACROSTRAINEXTENDED_H_
//! Substrate boundary condition  
class MacrostrainExtended : public MacrostrainBoundaryProperties
{
 public:
  
  virtual  ~MacrostrainExtended() {};
 

  static MacrostrainExtended* create(const ModelOptions& options);

 protected:

  MacrostrainExtended(const ModelOptions& options);

  virtual void 	do_init (void);

 private:

  
  

};


inline
MacrostrainExtended::MacrostrainExtended(const ModelOptions& options)
 : MacrostrainBoundaryProperties(options)
{
}


inline MacrostrainExtended* MacrostrainExtended::create(const ModelOptions& options)
{
  return new MacrostrainExtended(options);
}




#endif
