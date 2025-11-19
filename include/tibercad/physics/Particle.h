// $Id$


//! A class to represent a particle population
/*!
 * The particle population is represented by the 
 * type of particle (its charge), its density,
 * an equivalent quasi-Fermi level according to 
 * a thermal equilibrium Fermi-Dirac statistics,
 * an effective carrier temperature and by the
 * derivative with respect to the quasi Fermi level.
 *
 * This class is mainly thought for the exchange of
 * density data between models.
 */
class Particle
{

  public:

    //! Constructor
    Particle(double charge, double density, double Ef, double kT);

    //! Get the density
    double density(void) const
    { return _density; }

    //! Get the density
    double kT(void) const
    { return _kT; }

    //! Get the quasi Fermi level
    double fermi_level(void) const
    { return _fermi_level; }


  private:

    //! The charge in units of the elementary charge
    double _charge;

    //! The density in cm^-3
    double _density;

    //! The quasi Fermi level
    double _fermi_level;

    //! The effective temperature in eV
    double _kT;

};

inline
Particle::Particle(double charge, double density, double Ef, double kT) :
  _charge(charge),
  _density(density),
  _fermi_level(Ef),
  _kT(kT)
{
}
