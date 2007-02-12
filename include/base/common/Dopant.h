// $Id$

#ifndef _DOPANT_H_
#define _DOPANT_H_


//! Describes a dopant with a single energy level
class Dopant
{

  public:

    //! The type of the dopant
    enum DopingType
    {
      P_TYPE = -1, /*!< p-type, acceptor */
      N_TYPE = 1   /*!< n-type, donor */
    };

    //! Constructor
    /*!
     * \param density doping density
     * \param ionisation_energy the distance of the doping level from
     * the corresponding band edge \f$\vert E_{c,v} - E_d\vert\f$
     * \param g_factor g (cf. get_ionized_dopant_density())
     * \param type the doping type
     */
    Dopant(double density = 0.0, double ionisation_energy = 0.025,
        int g_factor = 2, DopingType type = N_TYPE);


    //! Copy constructor
    Dopant(const Dopant& dopant);

    
    //! Get the doping density
    double get_doping_density(void) const;

    
    //! Get the ionisation energy of the dopant
    /*!
     * \f$E_i = \vert E_{c,v} - E_d\vert\f$
     */
    double get_ionisation_energy(void) const;

    
    //! Get g (cf. get_ionized_dopant_density())
    int get_g_factor(void) const;

    
    //! Get the doping type
    DopingType get_type(void) const;

    
    void set_doping_density(double N);

    void set_ionisation_energy(double E_i);

    void set_g_factor(int g);

    void set_type(DopingType type);

    //! Get the density of ionized dopants
    /*!
     * For donors this is
     * \f[
     * N_d^+ = \frac{N_d}{1 + g \exp{??}}
     * \f]
     *
     */
    double get_ionized_dopant_density(double arg, double kT);


    /*! \brief {Get the derivative of the ionized doping density with
     * respect to the electric potential}
     */
    double get_ionized_dopant_density_derivative(double arg, double kT);


  private:

    //! The doping density
    double _density;

    //! The doping type
    DopingType _type;

    //! The ionisation energy
    double _ionisation_energy;

    //! The g factor
    int _g_factor;

};



//
// inline member functions
//

inline
Dopant::Dopant(double density, double ionisation_energy,
               int g_factor, DopingType type)
  : _density(density),
    _type(type),
    _ionisation_energy(ionisation_energy),
    _g_factor(g_factor)
{
}

inline
Dopant::Dopant(const Dopant& dopant)
  : _density(dopant._density),
    _type(dopant._type),
    _ionisation_energy(dopant._ionisation_energy),
    _g_factor(dopant._g_factor)
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
}

inline
void
Dopant::set_ionisation_energy(double E_i)
{
  _ionisation_energy = E_i;
}

inline
void
Dopant::set_g_factor(int g)
{
  _g_factor = g;
}

inline
void
Dopant::set_type(Dopant::DopingType type)
{
  _type = type;
}

#endif //_DOPANT_H_
