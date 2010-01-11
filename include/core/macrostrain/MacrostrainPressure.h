#ifndef _MACROSTRAINPRESSURE_H_
#define _MACROSTRAINPRESSURE_H_

//!mechanical pressure boundary condition class
class MacrostrainPressure : public MacrostrainBoundaryProperties
{
 public:


  double get_value(void) const;

  static MacrostrainPressure* create(const ModelOptions& options);

 protected:

  MacrostrainPressure(const ModelOptions& options);

  virtual void 	do_init (void);

  


 private:

  //! value of pressure in GPa
  double value;

}; 

inline MacrostrainPressure* MacrostrainPressure::create(const ModelOptions& options)
{
  return new MacrostrainPressure(options);
}


inline double MacrostrainPressure::get_value(void) const
{
  return value;
}
#endif
