// $Id$

#ifndef _PIEZOPOLARIZATION_H_
#define _PIEZOPOLARIZATION_H_

#include "tibercad/physics/misc/PolarizationModel.h"
#include "tibercad/module/SolutionProvider.h"
#include "tibercad/base/tiber_dll.h"


// Basic Piezopolarization model
class TBDLLOCAL Piezopolarization: public PolarizationModel
{

  public:

    virtual ~Piezopolarization(void) {};

    static Piezopolarization* create(const ModelOptions& options);


  protected:

    Piezopolarization(const ModelOptions& options);

    virtual void do_init(void);

    virtual void read_database(void);

    virtual void do_calculate(const libMesh::Elem* elem, const libMesh::Point& point);

  private:

    //! The strain simulation
    SolutionProvider _strain;

    //! Piezoelectric modulus \f$e_{33}\f$ (wurtzite)
    double _e33;

    //! Piezoelectric modulus \f$e_{31}\f$ (wurtzite)
    double _e31;

    //! Piezoelectric modulus \f$e_{15}\f$ (wurtzite) or \f$e_{14}\f$ (zincblende)
    union
    {
        double _e15;
        double _e14;
    };
};


inline
Piezopolarization::Piezopolarization(const ModelOptions& options) :
  PolarizationModel(options),
  _e33(0),
  _e31(0),
  _e15(0)
{
}


inline
Piezopolarization*
Piezopolarization::create(const ModelOptions& options)
{
  return new Piezopolarization(options);
}


#endif // _PIEZOPOLARIZATION_H_
