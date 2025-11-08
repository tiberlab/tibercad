#ifndef _ZB_FREEDYNAMICALMATRIX_H_
#define _ZB_FREEDYNAMICALMATRIX_H_


#include "PhononModel.h"
#include "DynamicalMatrix.h"

class ZbFreeDynamicalMatrix: public DynamicalMatrix
{
 public:
  //!constructor
  ZbFreeDynamicalMatrix(const ModelOptions& options) : DynamicalMatrix(options) {};

  //!destructor
  ~ZbFreeDynamicalMatrix() {};

  //! Create a ZbLatticeThermalConductivity object
 static  ZbFreeDynamicalMatrix* create(const ModelOptions& options);


//! Update the lattice thermal conductivity given the Temperature
  virtual void re_init(); 


  virtual void set_phonon_model(PhononModel* phonon_model);
 
 private:

 PhononModel* _phonon_model;
 double w0;

 protected:

  virtual void read_database(void);

  virtual void do_init(void);

  inline  virtual PhysicalModel*  create_new (void) const;

};

inline
ZbFreeDynamicalMatrix* ZbFreeDynamicalMatrix::create(const ModelOptions& options)
{
  return (new ZbFreeDynamicalMatrix(options));
}

inline
PhysicalModel*  ZbFreeDynamicalMatrix::create_new (void) const
{
  return (new  ZbFreeDynamicalMatrix(get_options()) );
}

inline
void ZbFreeDynamicalMatrix::set_phonon_model(PhononModel* phonon_model)
{
  _phonon_model = phonon_model;
}



#endif
