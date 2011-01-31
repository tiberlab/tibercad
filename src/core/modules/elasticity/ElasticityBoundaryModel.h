// $Id$

#ifndef _ELASTICITYBOUNDARYMODEL_H_tens_
#define _ELASTICITYBOUNDARYMODEL_H_tens_

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
		   const Point& point);

    //! Creator function
    static ElasticityBoundaryModel* create(const ModelOptions& options);

    const void get_coefficients(RealTensor& H, RealGradient& R);

    void set_normal(const Point p);

  protected:

    //! Constructor
    ElasticityBoundaryModel(const ModelOptions& options);

    virtual void do_init(void);

    virtual void create_submodels(void){};

  Point _normal;

  private:

  BoundaryModel* _bm;

  static TiberModelObject*  _create(const ModelOptions& options);
  
  static void  _destroy( TiberModelObject* p);

  RealTensor _H_tens;

  RealGradient _R_vec;

  double _A_scal;
   
};

inline
const
void
ElasticityBoundaryModel::get_coefficients(RealTensor& H, RealGradient& R)
{
  H = _H_tens;
  R = _R_vec;
}

inline
void
ElasticityBoundaryModel::set_normal(const Point normal)
{
  _normal = normal;
}



inline
ElasticityBoundaryModel::ElasticityBoundaryModel(const ModelOptions& options) :
  PhysicalModel(options),
  _H_tens(0),
  _R_vec(0),
  _A_scal(0)
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

#endif // _MYPOISSONMODEL_H_tens_
