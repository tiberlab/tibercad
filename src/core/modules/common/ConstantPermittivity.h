// $Id$

#ifndef _CONSTANTPERMITTIVITY_H_
#define _CONSTANTPERMITTIVITY_H_

#include "PermittivityModel.h"
#include "SimulationInterface.h"
#include "tiber_dll.h"
#include "Database.h"
class Elem;
class Point;

// Base class for charge density models
class  TBDLLOCAL ConstantPermittivity: public PermittivityModel
{

  public:

    virtual ~ConstantPermittivity(void) {};

    static ConstantPermittivity* create(const ModelOptions& options);


  protected:

    ConstantPermittivity(const ModelOptions& options);

    virtual void do_init(void);

    virtual void calculate(const Elem* elem, const Point& point){};

    virtual void read_database(void);

  private:

    RealVectorValue _permittivity_diag;

};


inline
ConstantPermittivity::ConstantPermittivity(const ModelOptions& options) :
  PermittivityModel(options),
  _permittivity_diag(0)
{
}


inline
ConstantPermittivity*
ConstantPermittivity::create(const ModelOptions& options)
{
  return new ConstantPermittivity(options);
}


#endif // _PIEZOPOLARIZATION_H_
