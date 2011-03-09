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

  const bool is_extended(void);

  const void get_coefficients(RealTensor& H, double& A,RealGradient& R);

    void set_normal(const Point p);


  protected:

    //! Constructor
    ElasticityBoundaryModel(const ModelOptions& options);

  void set_coefficients(RealTensor H,double A, RealGradient R);
  
  void set_is_extended(bool is_extended);

  Point _normal; 

  private:


 

  //! constrain matrix
  RealTensor _H_tens; 
 
  //! is extended
  double _coeff;

  //! contrain vector
  RealGradient _R_vec;

  //! is extended
  bool _is_extended;

  static TiberModelObject*  _create(const ModelOptions& options);
  
  static void  _destroy( TiberModelObject* p);
 
};

inline
const
void
ElasticityBoundaryModel::get_coefficients(RealTensor& H, double& A, RealGradient& R)
{
  H = _H_tens;
  R = _R_vec;
  A = _coeff;
}

inline
void
ElasticityBoundaryModel::set_normal(const Point normal)
{
  _normal = normal;
}


inline
void
ElasticityBoundaryModel::set_is_extended(bool is_extended)
{
  _is_extended = is_extended;
}

inline
const
bool 
ElasticityBoundaryModel::is_extended(void)
{
  return _is_extended;
}

inline
void
ElasticityBoundaryModel::set_coefficients(RealTensor H, double A, RealGradient R)
{
  _H_tens = H;
  _R_vec = R;
  _coeff = A;
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
