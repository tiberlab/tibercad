#ifndef _ZB_RAMANTENSOR_H_
#define _ZB_RAMANTENSOR_H_


#include "PhononModel.h"
#include "RamanTensor.h"
#include "Macrostrain.h"


class ZbRamanTensor: public RamanTensor
{
 public:
  //!constructor
  ZbRamanTensor(const ModelOptions& options) ;

  //!destructor
  ~ZbRamanTensor() {};

  //! Create a ZbLatticeThermalConductivity object
  static    ZbRamanTensor* create(const ModelOptions& options);

 
  virtual void re_init(){}; 

  virtual void set_phonon_model(PhononModel* phonon_model);
 
 private:

 PhononModel* _phonon_model;

 double _raman_d;
 
  protected:

  virtual void read_database(void);

  virtual void do_init(void);

  inline  virtual PhysicalModelInterface*  create_new (void) const;


};

inline
ZbRamanTensor* ZbRamanTensor::create(const ModelOptions& options)
{
  return (new ZbRamanTensor(options));
}

inline
PhysicalModelInterface*  ZbRamanTensor::create_new (void) const
{
  return (new  ZbRamanTensor(get_options()) );
}

inline
void  ZbRamanTensor::set_phonon_model(PhononModel* phonon_model)
{
  _phonon_model = phonon_model;
}



#endif
