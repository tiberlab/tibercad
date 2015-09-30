// $Id: DensityOfStates.h 2117 2010-10-20 20:56:13Z maufder $

#ifndef _DENSITYOFSTATES_H_
#define _DENSITYOFSTATES_H_


#include "PhysicalModelInterface.h"


/*!
 * \brief Base class for density of states
 */
class TBDLEXPORT DensityOfStates : public PhysicalModelInterface
{

  public:

    //! Destructor
    virtual ~DensityOfStates(void) {};


    //! Creator function
    static DensityOfStates* create(const ModelOptions& options);


    //! Set the reference energy
    /*!
     * \note some implementations may ignore the reference energy!
     */
    void set_reference_energy(double E0);



    //! Get occupied states
    /*!
     * \return the density of occupied states in cm^-3
     * for the given quasi Fermi-level \c Ef and the derivative
     *
     * \param Ef quasi Fermi-level in eV
     * \param Epot potential energy (electrostatic energy) in eV
     * \param kT temperature in eV
     */
    std::pair<double, double>
    get_occupied_density_and_derivative(double Ef, double Epot,
        double kT, const Elem* elem, const Point& p, double kTlattice = -1) const;

    //overloading for Trap.C
    std::pair<double, double>
    get_occupied_density_and_derivative(double Ef, double Epot,
        double kT, double kTlattice = -1) const;

    //! Get the particle name
    char get_particle(void) const;


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


    //! The reference energy
    const std::vector<double>& get_reference_energy(void) const
        { return _reference_energy; }

    //! The effective mass
    const std::vector<double>& get_effective_mass(void) const
        { return _effective_mass; }

    //! The reference energy
    std::vector<double>& reference_energy(void)
        { return _reference_energy; }

    //! The effective mass
    std::vector<double>& effective_mass(void)
        { return _effective_mass; }


    //! The effective DOS
    double get_effective_dos(void) const { return _effective_dos; }


  protected:

    //! Constructor
    DensityOfStates(const ModelOptions& options);

    //! A readable reference to the effective DOS
    double& effective_dos(void) { return _effective_dos; }

    //! Calculate density and derivative
    virtual std::pair<double, double>
    calculate_density_and_derivative(double Ef, double Epot, 
              double kT, double kTlattice, const Elem* elem, const Point& p) const = 0;

    //overloading for Trap.C
    virtual std::pair<double, double>
    calculate_density_and_derivative(double Ef, double Epot, double kT, double kTlattice) const = 0;

    // Do we have a quantum density?
    bool& is_quantum_density(void);



  private:

    //! The reference energy
    std::vector<double> _reference_energy;

    //! The effective_mass
    std::vector<double> _effective_mass;

    //! The effective DOS
    double _effective_dos;

    //! We keep the name of the particle
    /*!
     * e -> electron
     * h -> hole
     */
    char _particle;


    //! Do or don't use quantum density
    bool _use_quantum;


    //! \c true if the last calculated density is a quantum density
    bool _is_quantum;

};

//
// inline methods
//

inline
void
DensityOfStates::set_reference_energy(double E0)
{
  _reference_energy[0] = E0;
}


//inline
//void
//DensityOfStates::set_effective_dos(double Neff)
//{
//  _effective_dos = Neff;
//}

inline
char
DensityOfStates::get_particle(void) const
{
  return _particle;
}

inline
std::pair<double, double>
DensityOfStates::get_occupied_density_and_derivative(double Ef, double Epot,
    double kT, const Elem* elem, const Point& p, double kTlattice) const
{
  kTlattice = (kTlattice < 0) ? kT : kTlattice;
  return calculate_density_and_derivative(Ef, Epot, kT, kTlattice, elem, p);
}

inline
std::pair<double, double>
DensityOfStates::get_occupied_density_and_derivative(double Ef, double Epot,
    double kT, double kTlattice) const
{
  kTlattice = (kTlattice < 0) ? kT : kTlattice;
  return calculate_density_and_derivative(Ef, Epot, kT, kTlattice);
}

inline
void
DensityOfStates::use_quantum_density(bool use_quantum)
{
  _use_quantum = use_quantum;
}


inline
bool
DensityOfStates::is_quantum_density(void) const
{
  return (_use_quantum & _is_quantum);
}

inline
bool&
DensityOfStates::is_quantum_density(void)
{
  return _is_quantum;
}


inline
bool
DensityOfStates::has_quantum_density(void) const
{
  return _use_quantum;
}


#endif // _DENSITYOFSTATES_H_
