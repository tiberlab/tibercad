// $Id$

#ifndef _FOURIERMODEL_H_
#define _FOURIERMODEL_H_

#include "HeatTransportModel.h"

#include "point.h"
#include "tensor_value.h"
#include "vector_value.h"
#include "PhysicalModelInterface.h"
#include "tiber_dll.h"

class Elem;


//! The base class for Poisson boundary conditions
class TBDLLOCAL FourierModel : public HeatTransportModel
{

  public:
 
    //! Creator function
    static FourierModel* create(const ModelOptions& options);

   //! Destructor
    ~FourierModel(void) {};

  protected:

    //! Initialize
    virtual void do_init(void);

 
    /* In some cases it might be useful to reimplement this: */
    // virtual void do_init_interface(const PhysicalModelInterface* comp_A,
    //         const PhysicalModelInterface* comp_B);


    /* This is not used here: */
    // virtual void read_database(void);


    /* We do not use this here: */
    // virtual void read_interface_database(void);


    //! Create a new object of the same type
    virtual PhysicalModelInterface* create_new(void) const;


  private:

   //! Constructor
  FourierModel(const ModelOptions& options);
  
};


inline
PhysicalModelInterface*
FourierModel::create_new(void) const
{
  return new  FourierModel(get_options());
}

inline
FourierModel*
FourierModel::create(const ModelOptions& options)
{
  return new  FourierModel(options);
}
#endif // _GRAYMODEL_H_
