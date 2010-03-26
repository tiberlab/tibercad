// $Id$

#ifndef _MYPOISSONMODEL_H_
#define _MYPOISSONMODEL_H_

#include "PhysicalModel.h"

#include "point.h"
#include "tensor_value.h"
#include "vector_value.h"


#undef TIBER_MODULE_PREFIX
#define TIBER_MODULE_PREFIX mypoisson

class Elem;

//! This is the base class for the Poisson physical model
class MyPoissonModel : public PhysicalModel
{

  public:

    //! Destructor
    virtual ~MyPoissonModel(void) {};

    //! Creator function
    static MyPoissonModel* create(const ModelOptions& options);


    //! Get the relative permittivity
    const RealTensor& get_permittivity(void) const;


    //! Get the total polarization
    const RealVectorValue& get_polarization(void) const;


    //! Get the total charge
    double get_charge_density(void) const;


    //! Set the relative permittivity
    RealTensor& get_permittivity(void);


    //! Set the total polarization
    RealVectorValue& get_polarization(void);


    //! Set the total charge
    void set_charge_density(double charge_density);

    //! Reinit for element \c elem
    void set_element(const Elem* elem);


    //! Set the current point
    void set_point(const Point& point);


    //! Calculate everything
    void calculate(void);


  protected:

    //! Constructor
    MyPoissonModel(const ModelOptions& options) : PhysicalModel(options) {}


    //! do the actual calculation
    virtual void do_calculate(void) = 0;



  private:


    //! The element we are currently working on
    const Elem* _elem;


    //! The point we are currently using
    Point _point;


    //! The relative permittivity
    RealTensor _permittivity;

    //! The polarization
    RealVectorValue _polarization;

    //! The charge density
    double _charge;

};


inline
void
MyPoissonModel::set_element(const Elem* elem)
{
  _elem = elem;
}


inline
void
MyPoissonModel::set_point(const Point& point)
{
  _point = point;
}


inline
const RealTensor&
MyPoissonModel::get_permittivity(void) const
{
  return _permittivity;
}

inline
RealTensor&
MyPoissonModel::get_permittivity(void)
{
  return _permittivity;
}


inline
RealVectorValue&
MyPoissonModel::get_polarization(void)
{
  return _polarization;
}


inline
const RealVectorValue&
MyPoissonModel::get_polarization(void) const
{
  return _polarization;
}


inline
double
MyPoissonModel::get_charge_density(void) const
{
  return _charge;
}

inline
void
MyPoissonModel::set_charge_density(double charge_density)
{
  _charge = charge_density;
}


inline
void
MyPoissonModel::calculate(void)
{
  do_calculate();
}

#endif // _MYPOISSONMODEL_H_
