#ifndef _MACROSTRAINSUBSTRATE_H_
#define _MACROSTRAINSUBSTRATE_H_
//! Substrate boundary condition  
class MacrostrainSubstrate : public MacrostrainBoundaryProperties
{
 public:
  
  ~MacrostrainSubstrate();

  Material* get_material(void) const;

  static MacrostrainSubstrate* create(const ModelOptions& options);

 protected:

  MacrostrainSubstrate(const ModelOptions& options);

  virtual void 	do_init (void);


 private:

  Material*  material;
  

};


inline MacrostrainSubstrate* MacrostrainSubstrate::create(const ModelOptions& options)
{
  return new MacrostrainSubstrate(options);
}
//! returns a pointer to the substrate material 
inline  Material* MacrostrainSubstrate::get_material( ) const
{

  return material;

}

#endif
