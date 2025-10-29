// $Id$

#ifndef _CLAMP_H_
#define _CLAMP_H_

#include "ElasticityBoundaryModel.h"

//class Elem;


//! The base class for Poisson boundary conditions
class TBDLLOCAL Clamp : public ElasticityBoundaryModel
{

  public:

    //! Destructor
    ~Clamp(void) {};

    //! Creator function
    static Clamp* create(const ModelOptions& options);


    //! Calculate for a point on the given side
    virtual void calculate(const libMesh::Elem* elem, unsigned int side,
			   const libMesh::Point& point){};

  protected:

    //! Initialize
  virtual void do_init(void);


  private:

    //! Constructor
    Clamp(const ModelOptions& options);

};



inline
Clamp::Clamp(const ModelOptions& options) :
  ElasticityBoundaryModel(options)
{
}


inline
Clamp*
Clamp::create(const ModelOptions& options)
{
  return new Clamp(options);
}


#endif // _POISSONDIRICHLET_H_
