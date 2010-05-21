// $Id$

#ifndef _TRAP_H_
#define _TRAP_H_


#include "PhysicalModelInterface.h"


/*!
 * \brief Base class for traps in semiconductors
 *
 * The implemented traps in this model are:
 * \li neutral electron traps, type = ``eNeutral''
 * \li neutral hole traps, type = ``hNeutral''
 * \li donor like traps, type = ``donor''
 * \li acceptor like traps, type = ``acceptor''
 */
class Trap : public PhysicalModelInterface
{

  public:

    //! The trap type
    enum Type
    {
      NEUTRAL,  //!< neutral if unoccupied
      CHARGED   //!< charged if unoccupied
    };

    //! Destructor
    virtual ~Trap(void) {};


    //! Creator function
    static Trap* create(const ModelOptions& options);


    /*!
     * \brief Set the relevant energies
     *
     * \param Ec the conduction band edge in eV
     * \param Ev the valence band edge in eV
     * \param Ef the quasi fermi level in eV (\f$E_f = -q\phi\f$)
     * \param kT the thermal energy (in eV
     */
    void set_energies(double Ec, double Ev, double Ef, double kT);


    //! Get ionized density (= charge density, cm^-3)
    /*!
     * \return the charge density in cm^-3 including sign
     */
    double get_ionized_density(void) const;


    //! Get the derivative with respect to the quasi fermi level
    double get_ionized_density_derivative(void) const;


    //! Get the particle type
    char get_particle(void) const;


    //! Get the trap type
    Type get_type(void) const;


  protected:

    //! Constructor
    Trap(const ModelOptions& options);


    //! Initialize
    virtual void do_init(void);


    //! Set the trap type
    void set_type(Type type);


  private:

    //! The trap density in cm^-3
    double _density;

    //! The trap type
    Type _type;

    //! The particle type
    char _particle;

    //! The trap level
    double _level;

    //! The energy reference
    char _energy_reference;


    //! The conduction band
    double _Ec;

    //! The valence band
    double _Ev;

    //! The quasi fermi level
    double _fermi_level;

    //! The thermal energy
    double _kT;


    //! Calculate the trap level
    double _trap_level(void) const;

};

//
// inline methods
//


inline
Trap*
Trap::create(const ModelOptions& options)
{
  return new Trap(options);
}



inline
void
Trap::set_energies(double Ec, double Ev, double Ef, double kT)
{
  _Ec = Ec;
  _Ev = Ev;
  _fermi_level = Ef;
  _kT = kT;
}


inline
char
Trap::get_particle(void) const
{
  return _particle;
}


inline
Trap::Type
Trap::get_type(void) const
{
  return _type;
}


inline
void
Trap::set_type(Type type)
{
  _type = type;
}

#endif // _TRAP_H_
