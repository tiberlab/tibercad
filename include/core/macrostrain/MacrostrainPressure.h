#ifndef _MACROSTRAINPRESSURE_H_
#define _MACROSTRAINPRESSURE_H_

//!mechanical pressure boundary condition class
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
#endif
