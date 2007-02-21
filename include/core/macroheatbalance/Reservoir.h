#ifndef _RESERVOIR_H_
#define _RESERVOIR_H_

#include "ThermalContact.h"



class Reservoir : public ThermalContact
{
 public:
  Reservoir(){};
  
  ~Reservoir(){};

  double get_temperature(void) const;

  static Reservoir* create(void);

 protected:

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
