// $Id$

#ifndef _PARTICLEDENSITY_H_
#define _PARTICLEDENSITY_H_

#include "PhysicalModelInterface.h"
#include "TiberCad.h"
#include "TypeDefs.h"
#include "Constants.h"
#include "IDSet.h"
#include "SimulationOptions.h"
#include "tiber_dll.h"

#include "point.h"

#include <string>
#include <set>
#include <vector>

class SimulationInterface;
class Embracing;
class Elem;


//! An abstract class for a classical particle density
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
 * The length units are cm.
 *
 */
class ParticleDensity : public PhysicalModelInterface
{

  public:

    //! The minimum density
    static const double MINDENSITY;


    //! Constructor
    ParticleDensity(const ModelOptions& options);


    /*!
     * \param name the particle name (electron, hole etc.)
     * \param statistics the statistics to be used for classical densities
     */
    ParticleDensity(const std::string& name,
        TiberCad::Statistics statistics = TiberCad::BOLTZMANN);


    //! Create a ParticleDensity object
    static ParticleDensity* create(const ModelOptions& options);


    // Destructor
    ~ParticleDensity(void) {};


    //! Get the particle name
    const std::string& get_particle_name(void);


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
     * \param mDOS the DOS effective mass
     * \param E the particle energy
     * \param E_F the electro-chemical potential
     * \param kT the temperature in eV
     */
    void set_classical_parameters(double mDOS, double E,
        double E_F = 0.0, double kT = -1.0);


    //! Set the parameters for the classical density
    /*!
     * For compatibility only.
     *
     * \param Neff the effective DOS
     * \param E the particle energy
     * \param E_F the electro-chemical potential
     * \param kT the temperature in eV
     */
    void set_classical_parameters_comp(double Neff, double E,
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
     * \return the particle density in cm^-3
     *
     * \attention 
     *   The particle density is calculated only once for given parameters.
     *   A call to a method that changes parameters will force a
     *   recalculation at the next call.
     * 
     */
    double get_particle_density(void);


    //! Get particle density derivative
    /*!
     * set_classical_parameters() and set_element_and_point() have to be
     * called before this one. The derivative is given with respect to the
     * electrostatic potential/quasi Fermi level.
     *
     * \attention 
     *   The particle density derivative is calculated only once for given
     *   parameters. A call to a method that changes parameters will force a
     *   recalculation at the next call.
     * 
     */
    double get_particle_density_derivative(void);


    //! Get \f$\gamma = F_{1/2}(\mu)/exp(\mu)\f$
    double get_gamma(void);


    //! Set up an embracing region
    void set_embracing(Embracing* embracing);


    //! Get a pointer to the quantum density simulation
    SimulationInterface* get_quantum_simulation(void);



  protected:

    //! Do initialization
    virtual void do_init(void);



  private:

    //! The name of the particle
    std::string _name;

    //! The charge per particle in units of the electron charge
    double _charge;

    
    //! The statistics
    TiberCad::Statistics _statistics;


    //! Do or don't use quantum density
    bool _use_quantum;


    //! \c true if the last calculated density is a quantum density
    bool _is_quantum;


    //! The quantum density calculation
    std::vector<SimulationInterface*> _quantum_density;


    //! The region IDs of the barrier
    IDSet _barrier_ids;


    //! The ID of the Density variable
    std::vector<ID> _density_ids;


    //! The ID for the BandEdge3D variable
    std::vector<ID> _3D_edge;


    //! The element we are currently working on
    const Elem* _elem;


    //! The point we are currently working on
    Point _p;


    //! The effective density of states
    double _N_eff;

    //! The DOS mass, including degeneracy
    double _dos_mass;


    //! The DOS factor
    double _dos_factor;


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


    //! The parameter \f$\gamma = F_{1/2}(\mu)/exp(\mu)\f$
    double _gamma;


    //! The embracing of classical and quantum calculation
    Embracing* _embracing;


    //! If \c true a continuum will be added on top of the quantum density
    bool _add_continuum;


    //! Calculate the particle density
    void calculate_density(void);


    //! Calculate the quantum density
    /*!
     * returns \true if the quantum density could be evaluated
     */
    bool quantum_density(void) TBDLLOCAL;


    //! Calculate classical particle density
    template <TiberCad::Statistics>
    void classical_density(void) TBDLLOCAL;


    //! Calculate the effective DOS from the DOS mass
    void calculate_effective_DOS(void);

};



//
// inline methods
//


inline
const std::string&
ParticleDensity::get_particle_name(void)
{
  return _name;
}


inline
void
ParticleDensity::set_particle_charge(double particle_charge)
{
  _charge = particle_charge;
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
  return _charge;
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
  _density = -1.0;
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
ParticleDensity::calculate_effective_DOS(void)
{
  _N_eff = _dos_factor * std::pow(_kT * _dos_mass, 1.5);
}




inline
void
ParticleDensity::set_classical_parameters(double mDOS, double E,
    double E_F, double kT)
{
  //_N_eff = N_eff;
  _dos_mass = mDOS;
  _E = E;
  _E_F = E_F;
  _kT = (kT >= 0) ? kT : SimulationOptions::temperature * Constants::k_B;

  calculate_effective_DOS();

  _argument = (_E_F - E) / _kT;

  _density = -1.0;
}


inline
void
ParticleDensity::set_classical_parameters_comp(double N_eff, double E,
    double E_F, double kT)
{
  _N_eff = N_eff;
  _dos_mass = 1; // TODO calculate DOS mass from N_eff
  _E = E;
  _E_F = E_F;
  _kT = (kT >= 0) ? kT : SimulationOptions::temperature * Constants::k_B;

  //calculate_effective_DOS();

  _argument = (_E_F - E) / _kT;

  _density = -1.0;
}






inline
void
ParticleDensity::set_element_and_point(const Elem* elem, const Point& p)
{
  _elem = elem;
  _p = p;

  _density = -1.0;
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
  if (_density < 0)
    calculate_density();

  return _density_derivative;
}



inline
double
ParticleDensity::get_gamma(void)
{
  if (_density < 0)
    calculate_density();

  return _gamma;
}



inline
SimulationInterface*
ParticleDensity::get_quantum_simulation(void)
{
  if (_quantum_density.size() > 0)
    return _quantum_density[0];

  return NULL;
}



#endif // _PARTICLEDENSITY_H_
