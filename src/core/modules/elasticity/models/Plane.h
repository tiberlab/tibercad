// $Id: Clamp.h 2124 2010-10-22 14:00:17Z gromano $

#ifndef _Plane_H_
#define _Plane_H_

#include "ElasticityBoundaryModel.h"

#include "point.h"
#include "tensor_value.h"
#include "vector_value.h"



class Elem;


//! The base class for Poisson boundary conditions
class TBDLLOCAL Plane : public ElasticityBoundaryModel
{

  public:

    //! Destructor
    ~Plane(void) {};

    //! Creator function
    static Plane* create(const ModelOptions& options);


    //! Calculate for a point on the given side
    virtual void calculate(const libMesh::Elem* elem, unsigned int side,
			   const libMesh::Point& point);


  protected:

    //! Initialize
  virtual void do_init(void){};

    /* In some cases it might be useful to reimplement this: */
    // virtual void do_init_interface(const PhysicalModelInterface* comp_A,
    //         const PhysicalModelInterface* comp_B);


    /* This is not used here: */
    // virtual void read_database(void);


    /* We do not use this here: */
    // virtual void read_interface_database(void);


  private:

    //! Constructor
    Plane(const ModelOptions& options);

   
};



inline
Plane::Plane(const ModelOptions& options) :
  ElasticityBoundaryModel(options)
{
}



inline
Plane*
Plane::create(const ModelOptions& options)
{
  return new Plane(options);
}




#endif // _POISSONDIRICHLET_H_
