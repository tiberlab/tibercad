// $Id$

#ifndef _ELASTICITYBOUNDARYMODEL_H_tens_
#define _ELASTICITYBOUNDARYMODEL_H_tens_

#include "PhysicalModel.h"

#include "point.h"
#include "tensor_value.h"
#include "vector_value.h"
#include "tiber_dll.h"

class Elem;

//! This is the base class for the Poisson physical model
class ElasticityBoundaryModel : public PhysicalModel
{

  public:

    //! Destructor
    virtual ~ElasticityBoundaryModel(void) {};
  
    //! Calculate properties 
    virtual void calculate(const Elem* elem, unsigned int side,
		   const Point& point) = 0;

    //! Creator function
    static ElasticityBoundaryModel* create(const ModelOptions& options);

    const void get_coefficients(RealTensor& H, RealGradient& R, double& A);

    void set_normal(const Point p);


  protected:

    //! Constructor
    ElasticityBoundaryModel(const ModelOptions& options);

    //virtual void do_init(void);

    //virtual void create_submodels(void){};

  void set_coefficients(RealTensor H, RealGradient R, double A);

  Point _normal; 

  private:


  double _is_extended;

  RealTensor _H_tens;

  RealGradient _R_vec;

  static TiberModelObject*  _create(const ModelOptions& options);
  
  static void  _destroy( TiberModelObject* p);
 
};

inline
const
void
ElasticityBoundaryModel::get_coefficients(RealTensor& H, RealGradient& R, double& A)
{
  H = _H_tens;
  R = _R_vec;
  A = _is_extended;
}

inline
void
ElasticityBoundaryModel::set_normal(const Point normal)
{
  _normal = normal;
}



inline
void
ElasticityBoundaryModel::set_coefficients(RealTensor H, RealGradient R, double A)
{
  _H_tens = H;
  _R_vec = R;
  _is_extended = A;
}

inline
ElasticityBoundaryModel::ElasticityBoundaryModel(const ModelOptions& options) :
  PhysicalModel(options),
  _H_tens(0),
  _R_vec(0),
  _is_extended(0)
{
}




//inline
//ElasticityBoundaryModel*
//ElasticityBoundaryModel::create(const ModelOptions& options)
//{
 
//  return dynamic_cast<ElasticityBoundaryModel*>(PhysicalModelInterface::create(_create,_destroy,options));
  
//}

#endif // _MYPOISSONMODEL_H_tens_
