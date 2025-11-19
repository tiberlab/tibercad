// $Id$

#ifndef _PYROPOLARIZATION_H_
#define _PYROPOLARIZATION_H_

#include "tibercad/model_base/PolarizationModel.h"
#include "tibercad/module/SimulationInterface.h"
#include "tibercad/base/tiber_dll.h"

// Base class for charge density models
class  TBDLLOCAL Pyropolarization: public PolarizationModel
{

  public:
  
   virtual ~Pyropolarization(void) {};
   
   static Pyropolarization* create(const ModelOptions& options);
  
  protected:

   Pyropolarization(const ModelOptions& options);

   virtual void do_init(void);

   virtual void read_database(void);

   virtual void do_calculate(const libMesh::Elem*, const libMesh::Point&) {};

  private:

    double _Pz;

    //! Initialize P from given Pz
    void _initP(void);

};


inline
Pyropolarization::Pyropolarization(const ModelOptions& options) :
  PolarizationModel(options),
  _Pz(0)
{
}


inline
Pyropolarization*
Pyropolarization::create(const ModelOptions& options)
{
  return new Pyropolarization(options);
}


#endif // _PIEZOPOLARIZATION_H_
