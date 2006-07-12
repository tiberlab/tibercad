// $Id$

#ifndef _DOPANT_H_
#define _DOPANT_H_

class Dopant
{

  public:

    enum DopingType
    {
      P_TYPE = -1,
      N_TYPE = 1
    };

    Dopant(double density = 0.0, double ionisation_energy = 0.025,
        int g_factor = 2, DopingType type = N_TYPE);

    double get_doping_density(void) const;

    double get_ionisation_energy(void) const;

    int get_g_factor(void) const;

    DopingType get_type(void) const;

    void set_doping_density(double N);

    void set_ionisation_energy(double E_i);

    void set_g_factor(int g);

    void set_type(DopingType type);

    double get_ionized_dopant_density(double arg, double kT);

    double get_ionized_dopant_density_derivative(double arg, double kT);


  private:

    DopingType _type;
    double _density;
    double _ionisation_energy;
    int _g_factor;

};


//
// inline member functions
//

inline
Dopant::Dopant(double density, double ionisation_energy,
               int g_factor, DopingType type)
  : _density(density), _type(type),
    _ionisation_energy(ionisation_energy), _g_factor(g_factor)
{
}

inline
double
Dopant::get_doping_density(void) const
{
  return _density;
}

inline
double
Dopant::get_ionisation_energy(void) const
{
  return _ionisation_energy;
}

inline
int
Dopant::get_g_factor(void) const
{
  return _g_factor;
}

inline
Dopant::DopingType
Dopant::get_type(void) const
{
  return _type;
}

inline
void
Dopant::set_doping_density(double N)
{
  _density = N;
};

inline
void
Dopant::set_ionisation_energy(double E_i)
{
  _ionisation_energy = E_i;
};

inline
void
Dopant::set_g_factor(int g)
{
  _g_factor = g;
};

inline
void
Dopant::set_type(Dopant::DopingType type)
{
  _type = type;
};

#endif //_DOPANT_H_
