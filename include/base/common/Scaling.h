// $Id$

#ifndef _SCALING_H_
#define _SCALING_H_

#include <cassert>

#include "Constants.h"

//! A class to define scaling parameters
/*!
 * This class holds the parameters for the scaling from physical units
 * to calculation units.
 */
class Scaling
{

  public:

    //! Default constructor
    /*!
     * Sets all scaling parameters to 1
     */
    Scaling(void);

    //! Copy constructor
    Scaling(const Scaling& scaling);

    //! Assignment operator
    Scaling& operator=(const Scaling& other);


    //! The type of scaling
    enum ScalingType
    {
      NONE,   /*!< all scaling parameters are 1 */
      UNITS,  /*!< units scaling, for drift-diffusion */
      DEMARI, /*!< De Mari scaling, for drift-diffusion */
      CUSTOM, /*!< other */
    };

    //! Get the scaling type
    ScalingType get_scaling_type(void) const;


    //! Get the scaling parameter for potentials / energies
    double get_potential_scaling(void) const;

    
    //! Get the scaling parameter for length units
    double get_length_scaling(void) const;


    //! Get the mesh units for calculation
    double get_calc_mesh_units(void) const;

    
    //! Get the scaling parameter for time
    double get_time_scaling(void) const;

    
    //! Get the scaling parameter for mobilities
    double get_mobility_scaling(void) const;


    //! Get the scaling parameter for densities
    double get_density_scaling(void) const;

    double get_recombination_scaling(void) const;

    double get_current_density_scaling(void) const;

    double get_lambda_squared(void) const;

    
    void set_scaling_type(ScalingType type);

    void set_potential_scaling(double phi0);

    void set_length_scaling(double x0);

    void set_calc_mesh_units(double units);
  
    void set_mobility_scaling(double mu0);

    void set_density_scaling(double C0);


  private:

    ScalingType _type;
    
    double _potential;
    double _length;
    double _mesh_units;
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
Scaling::get_calc_mesh_units(void) const
{
  return _mesh_units;
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
Scaling::set_calc_mesh_units(double units)
{
  assert(units != 0.0);

  _mesh_units = units;
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

inline
Scaling&
Scaling::operator=(const Scaling& other)
{
  _type = other._type;
  _potential = other._potential;
  _length = other._length;
  _mesh_units = other._mesh_units;
  _mobility = other._mobility;
  _density = other._density;
}


#endif //_SCALING_H_
