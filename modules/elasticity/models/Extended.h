// $Id: Clamp.h 2352 2011-02-18 12:55:55Z maufder $

#ifndef _EXTENDED_H_
#define _EXTENDED_H_

#include "ElasticityBoundaryModel.h"


class Elem;


//! The base class for Poisson boundary conditions
class TBDLLOCAL Extended : public ElasticityBoundaryModel
{

  public:

    //! Destructor
    ~Extended(void) {};

    //! Creator function
    static Extended* create(const ModelOptions& options);


    //! Calculate for a point on the given side
    virtual void calculate(const libMesh::Elem* elem, unsigned int side,
			   const libMesh::Point& point);

  protected:

    //! Initialize
  virtual void do_init(void);


  private:

    //! Constructor
    Extended(const ModelOptions& options);

};



inline
Extended::Extended(const ModelOptions& options) :
  ElasticityBoundaryModel(options)
{
}


inline
Extended*
Extended::create(const ModelOptions& options)
{
  return new Extended(options);
}


#endif // _POISSONDIRICHLET_H_
