// $Id$

#ifndef _CONVERSEPIEZO_H_
#define _CONVERSEPIEZO_H_

#include "BodyForceModel.h"

#include "tibercad/module/SimulationInterface.h"
#include "point.h"
#include "tensor_value.h"
#include "vector_value.h"
#include "tibercad/base/tiber_dll.h"

class Elem;



//! The base class for Poisson boundary conditions
class TBDLLOCAL ConversePiezo : public BodyForceModel
{

  public:
 
  //! Destructor
  ~ConversePiezo(void){};
  
  //! Creator function
  static ConversePiezo* create(const ModelOptions& options);
  
  virtual void calculate(const libMesh::Elem* elem, const libMesh::Point& point);

  protected:

    //! Initialize
    virtual void do_init(void);

    /* In some cases it might be useful to reimplement this: */
    // virtual void do_init_interface(const PhysicalModel* comp_A,
    //         const PhysicalModel* comp_B);


    /* This is not used here: */
     virtual void read_database(void);


    /* We do not use this here: */
    // virtual void read_interface_database(void);



  private:
  
  SimulationInterface* _simul;

  double _e33;
  double _e31;
  double _e15;
  
  ID ElFieldID;


  //! Constructor
  ConversePiezo(const ModelOptions& options);
  
};




inline
ConversePiezo*
ConversePiezo::create(const ModelOptions& options)
{ 
  return new  ConversePiezo(options);
}




#endif // _GRAYMODEL_H_
