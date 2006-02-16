// $Id$

#ifndef _SCALING_H_
#define _SCALING_H_

#include <cassert>

#include "Constants.h"

class Scaling
{

  public:

    Scaling(void);
    Scaling(const Scaling& scaling);

    enum ScalingType
    {
      NONE,
      UNITS,
      DEMARI,
      CUSTOM
    };

    ScalingType get_scaling_type(void) const;


    double get_potential_scaling(void) const;

    double get_length_scaling(void) const;

    double get_mobility_scaling(void) const;

    double get_density_scaling(void) const;

    double get_time_scaling(void) const;

    double get_recombination_scaling(void) const;

    double get_current_density_scaling(void) const;

    double get_lambda_squared(void) const;

    
    void set_scaling_type(ScalingType type);

    void set_potential_scaling(double phi0);

    void set_length_scaling(double x0);

    void set_mobility_scaling(double mu0);

    void set_density_scaling(double C0);


  private:

    ScalingType _type;
    
    double _potential;
    double _length;
    double _mobility;
    double _density;

};


//
// inline member functions
// 

inline
Scaling::ScalingType
Scaling::get_scaling_type(void) const
{
  return _type;
}

inline
double
Scaling::get_potential_scaling(void) const
{
  return _potential;
}

inline
double
Scaling::get_length_scaling(void) const
{
  return _length;
}

inline
double
Scaling::get_mobility_scaling(void) const
{
  return _mobility;
}

inline
double
Scaling::get_density_scaling(void) const
{
  return _density;
}

inline
double
Scaling::get_time_scaling(void) const
{
  return _length * _length / (_potential * _mobility);
}

inline
double
Scaling::get_recombination_scaling(void) const
{
  return _potential * _mobility * _density / (_length * _length);
}

inline
double
Scaling::get_current_density_scaling(void) const
{
  return _potential * _mobility * _density / _length;
}

inline
double
Scaling::get_lambda_squared(void) const
{
  double tmp = _length * _length * _density * Constants::e;
  return _potential / tmp;
}

inline
void
Scaling::set_scaling_type(Scaling::ScalingType type)
{
  _type = type;
}

inline
void
Scaling::set_potential_scaling(double phi0)
{
  assert(phi0 != 0.0);

  _potential = phi0;
}

inline
void
Scaling::set_length_scaling(double x0)
{
  assert(x0 != 0.0);

  _length = x0;
}

inline
void
Scaling::set_mobility_scaling(double mu0)
{
  assert(mu0 != 0.0);

  _mobility = mu0;
}

inline
void
Scaling::set_density_scaling(double C0)
{
  assert(C0 != 0.0);

  _density = C0;
}


#endif //_SCALING_H_
