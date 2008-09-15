// $Id$

#ifndef _PARTICLEDENSITY_H_
#define _PARTICLEDENSITY_H_

#include "TiberCad.h"
#include "TypeDefs.h"
#include "Constants.h"
#include "SimulationOptions.h"

#include "point.h"

#include <string>
#include <set>

class SimulationInterface;
class Embracing;
class Elem;


//! An abstract class for a particle density
/*!
 *
 * The classical density is calculated as
 * \f[
 *   n = N_{eff}f(E_F - E)
 * \f]
 * where \f$N_{eff}\f$ is the effective density of states, \f$E_F\f$ the
 * electro-chemical potential and \f$E\f$ the particle energy.
 * The function \f$f\f$ can be the exponential or the Fermi function of order
 * 1/2, according to carrier statistics.
 *
 */
class ParticleDensity
{

  public:

    //! Constructor
    /*!
     * \param particle_charge the charge per particle, in units of the
     * electron charge \it e
     * \param statistics the statistics to be used for classical densities
     */
    ParticleDensity(double particle_charge = -1.0,
        TiberCad::Statistics statistics = TiberCad::BOLTZMANN);


    // Destructor
    ~ParticleDensity(void) {};


    //! Set the charge per particle
    /*!
     * The charge per particle has to be given in units of the electron
     * charge \it e, so e.g. for an electron -1, for a hole +1 etc.
     */
    void set_particle_charge(double particle_charge);


    //! Set the particle statistics
    void set_statistics(TiberCad::Statistics statistics);


    //! Add the name of a quantum density calculation for this particle
    /*!
     * Instead of calculating the statistical classical densities one can
     * use a quantum density calculation. There can be more than one 
     * quantum density calculation for the same particle as one could consider
     * eg. different valleys using different models.
     *
     * \param name the name of the quantum density calculation
     */
    void add_quantum_density(const std::string& name);


    //! Do or don't use the quantum density
    /*!
     * Independently on this flag classical densities will be calculated
     * when either there is no quantum density calculation or it wasn't
     * solved yet.
     */
    void use_quantum_density(bool use_quantum = true);


    // Do we have a quantum density?
    bool is_quantum_density(void) const;


    // Do we want to use a quantum density?
    bool has_quantum_density(void) const;


    //! Set the parameters for the classical density
    /*!
     * \param N_eff the effective density of states
     * \param E the particle energy
     * \param E_F the electro-chemical potential
     * \param kT the temperature in eV
     */
    void set_classical_parameters(double N_eff, double E,
        double E_F = 0.0, double kT = -1.0);


    //! Set the element and point we are working on
    void set_element_and_point(const Elem* elem, const Point& p);


    //! Get the particle charge
    double get_particle_charge(void) const;


    //! Get the statistics
    TiberCad::Statistics get_statistics(void) const;


    //! Get particle density
    /*!
     * set_classical_parameters() and set_element_and_point() have to be
     * called before this one.
     *
     * \attention {
     *   The particle density is calculated only once for given parameters.
     *   A call to a method that changes parameters will force a
     *   recalculation at the next call.
     * }
     */
    double get_particle_density(void);


    //! Get particle density derivative
    /*!
     * set_classical_parameters() and set_element_and_point() have to be
     * called before this one
     *
     * \attention {
     *   The particle density derivative is calculated only once for given
     *   parameters. A call to a method that changes parameters will force a
     *   recalculation at the next call.
     * }
     */
    double get_particle_density_derivative(void);


    //! Set up an embracing region
    void set_embracing(Embracing* embracing);


    //! Get a pointer to the quantum density simulation
    SimulationInterface* get_quantum_simulation(void);



  private:

    //! The charge per particle in units of the electron charge
    double _particle_charge;

    
    //! The statistics
    TiberCad::Statistics _statistics;


    //! Do or don't use quantum density
    bool _use_quantum;


    //! \c true if the last calculated density is a quantum density
    bool _is_quantum;


    //! The quantum density calculation
    SimulationInterface* _quantum_density;


    //! The ID of the density variable
    ID _density_id;


    //! The element we are currently working on
    const Elem* _elem;


    //! The point we are currently working on
    Point _p;


    //! The effective density of states
    double _N_eff;


    //! The electro-chemical potential in eV
    double _E_F;


    //! The particle energy in eV
    double _E;


    //! The temperature in eV
    double _kT;


    //! The argument that enters in the classical model
    /*!
     * \f$ \mathtt{_argument} = \frac{E_F - E}{kT}\f$
     */
    double _argument;


    //! The particle density
    double _density;


    //! The particle density derivative
    double _density_derivative;


    //! The embracing of classical and quantum calculation
    Embracing* _embracing;


    //! Calculate the particle density
    void calculate_density(void);


    //! Calculate the particle density
    /*!
     * Always call this after calculate_density()
     */
    void calculate_density_derivative(void);


    //! Calculate the quantum density
    /*!
     * returns \true if the quantum density could be evaluated
     */
    bool quantum_density(void);


    //! Calculate the quantum density derivative
    /*!
     * returns \true if the quantum density derivative could be evaluated
     *
     * \note For now, quantum density derivative is zero
     */
    bool quantum_density_derivative(void);


    //! Calculate classical particle density
    template <TiberCad::Statistics>
    void classical_density(void);


    //! Calculate classical particle density
    template <TiberCad::Statistics>
    void classical_density_derivative(void);

};



//
// inline methods
//



inline
void
ParticleDensity::set_particle_charge(double particle_charge)
{
  _particle_charge = particle_charge;
}



inline
void
ParticleDensity::set_statistics(TiberCad::Statistics statistics)
{
  _statistics = statistics;
}



inline
double
ParticleDensity::get_particle_charge(void) const
{
  return _particle_charge;
}



inline
TiberCad::Statistics
ParticleDensity::get_statistics(void) const
{
  return _statistics;
}



inline
void
ParticleDensity::use_quantum_density(bool use_quantum)
{
  _use_quantum = use_quantum;
  _density = _density_derivative = -1.0;
}


inline
bool
ParticleDensity::is_quantum_density(void) const
{
  return (_use_quantum & _is_quantum);
}


inline
bool
ParticleDensity::has_quantum_density(void) const
{
  return _use_quantum;
}




inline
void
ParticleDensity::set_classical_parameters(double N_eff, double E,
    double E_F, double kT)
{
  _N_eff = N_eff;
  _E = E;
  _E_F = E_F;
  _kT = (kT >= 0) ? kT : SimulationOptions::temperature * Constants::k_B;

  _argument = (_E_F - E) / _kT;

  _density = _density_derivative = -1.0;
}



inline
void
ParticleDensity::set_element_and_point(const Elem* elem, const Point& p)
{
  _elem = elem;
  _p = p;

  _density = _density_derivative = -1.0;
}



inline
double
ParticleDensity::get_particle_density(void)
{
  if (_density < 0)
    calculate_density();

  return _density;
}



inline
double
ParticleDensity::get_particle_density_derivative(void)
{
  if (_density_derivative < 0)
    calculate_density_derivative();

  return _density_derivative;
}



inline
SimulationInterface*
ParticleDensity::get_quantum_simulation(void)
{
  return _quantum_density;
}



#endif // _PARTICLEDENSITY_H_
