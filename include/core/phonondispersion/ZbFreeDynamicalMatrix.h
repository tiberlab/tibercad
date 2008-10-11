#ifndef _ZB_FREEDYNAMICALMATRIX_H_
#define _ZB_FREEDYNAMICALMATRIX_H_


#include "PhononModel.h"
#include "DynamicalMatrix.h"

class ZbFreeDynamicalMatrix: public DynamicalMatrix
{
 public:
  //!constructor
  ZbFreeDynamicalMatrix() {};

  //!destructor
  ~ZbFreeDynamicalMatrix() {};

  //! Create a ZbLatticeThermalConductivity object
 static  ZbFreeDynamicalMatrix* create();


//! Update the lattice thermal conductivity given the Temperature
  virtual void re_init(); 


  virtual void set_phonon_model(PhononModel* phonon_model);
 
 private:

 PhononModel* _phonon_model;
 double w0;

 protected:

  virtual void read_database(void);

  virtual void do_init(void);

  inline  virtual PhysicalModelInterface*  create_new (void) const;

};

inline
ZbFreeDynamicalMatrix* ZbFreeDynamicalMatrix::create()
{
  return (new ZbFreeDynamicalMatrix());
}

inline
PhysicalModelInterface*  ZbFreeDynamicalMatrix::create_new (void) const
{
  return (new  ZbFreeDynamicalMatrix() ); 
}

inline
void ZbFreeDynamicalMatrix::set_phonon_model(PhononModel* phonon_model)
{
  _phonon_model = phonon_model;
}



#endif
