// $Id$

#ifndef _ELASTICITYBOUNDARYMODEL_H_
#define _ELASTICITYBOUNDARYMODEL_H_

#include "PhysicalModel.h"

#include "point.h"
#include "tensor_value.h"
#include "vector_value.h"
#include "BoundaryModel.h"
#include "tiber_dll.h"

class Elem;

//! This is the base class for the Poisson physical model
class TBDLLOCAL ElasticityBoundaryModel : public PhysicalModel
{

  public:

    //! Destructor
    virtual ~ElasticityBoundaryModel(void) {};
  
    //! Calculate properties 
    void calculate(const Elem* elem, unsigned int side,
		   const Point& point){};

    //! Creator function
    static ElasticityBoundaryModel* create(const ModelOptions& options);

    const void get_coefficients(RealTensor& H, double& A, RealGradient& R);
 
  protected:

    //! Constructor
    ElasticityBoundaryModel(const ModelOptions& options);

    virtual void do_init(void);

    virtual void create_submodels(void){};

  private:

  static TiberModelObject*  _create(const ModelOptions& options);
  
  static void  _destroy( TiberModelObject* p);

  RealTensor _H;

  RealGradient _R;

  double _A;
   
};

inline
const
void
ElasticityBoundaryModel::get_coefficients(RealTensor& H, double& A, RealGradient& R)
{
  H = _H;
  R = _R;
  A = _A;
}




inline
ElasticityBoundaryModel::ElasticityBoundaryModel(const ModelOptions& options) :
  PhysicalModel(options),
  _H(0),
  _R(0),
  _A(0)
{
}

inline
TiberModelObject*  ElasticityBoundaryModel::_create(const ModelOptions& options)
{

  return new ElasticityBoundaryModel(options);

}

inline
void  ElasticityBoundaryModel::_destroy( TiberModelObject* p)
{

  delete p;

}


inline
ElasticityBoundaryModel*
ElasticityBoundaryModel::create(const ModelOptions& options)
{
 
  return dynamic_cast<ElasticityBoundaryModel*>(PhysicalModelInterface::create(_create,_destroy,options));
  
}

#endif // _MYPOISSONMODEL_H_
