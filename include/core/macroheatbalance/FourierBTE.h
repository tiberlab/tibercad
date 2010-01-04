#ifndef _FOURIERBTE_H_
#define _FOURIERBTE_H_

#include "ThermalContact.h"
#include "Variable.h"
#include "TemperatureInterface.h"

class FourierBTE : public ThermalContact
{
 public:
  //!Constructor	 
  FourierBTE();
  //!Destructor 
  ~FourierBTE(){};
  //!Return the temperature of the contact
  double get_temperature(const Elem* elem, const Point&);

  //!Set the temperature of contact
  void set_temperature(double Temp){};
   
  //!Create a Reservoir object and return its pointer
  static  FourierBTE* create(void);

  //! The interface to the lattice temperature simulation
  TemperatureInterface _lattice_temp;

 //! The interface to the lattice temperature simulation
  TemperatureInterface _global_lattice_temp;

  bool _global_simulation;

 protected:
  //!Initialize the model
  virtual void 	do_init (void);

  /*! \copydoc Variable::set_variable_value() */
  virtual void set_variable_value(double value, ID id = 0){};
  
  
  /*! \copydoc Variable::get_variable_value() */
  virtual double get_variable_value(ID id = 0){};


 private:

  double _temperature; 

};


inline
FourierBTE* 
FourierBTE::create()
{
  return new  FourierBTE();
}


//inline 
//void  
//FourierBTE::set_temperature(double Temp ) 
//{

//   _temperature = Temp;

//}


#endif
