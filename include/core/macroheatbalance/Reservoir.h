#ifndef _RESERVOIR_H_
#define _RESERVOIR_H_

#include "ThermalContact.h"



class Reservoir : public ThermalContact
{
 public:
  //!Constructor	 
  Reservoir(){};
  //!Destructor 
  ~Reservoir(){};
  //!Return the temperature
  double get_temperature(void) const;
  //!Create a Reservoir object and return its pointer
  static Reservoir* create(void);

 protected:
  //!Initialize the model
  virtual void 	do_init (void);


 private:

  double _temperature; 

};


inline Reservoir* Reservoir::create()
{
  return new Reservoir();
}

inline double  Reservoir::get_temperature( ) const
{

  return _temperature;

}

#endif
