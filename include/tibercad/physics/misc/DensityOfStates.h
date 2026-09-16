/*
 * This file is part of tiberCAD.
 *
 * tiberCAD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * tiberCAD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with tiberCAD. If not, see <https://www.gnu.org/licenses/>.
 */

/*!
 * \file DensityOfStates.h
 * \brief Public tiberCAD API header.
 *
 * \note This header is part of the public tiberCAD API.
 *       API version: 3.5
 */


#ifndef TC_DENSITYOFSTATES_H
#define TC_DENSITYOFSTATES_H


#include "tibercad/physics/PhysicalModel.h"

class ExternalProfile;

/*!
 * \brief Base class for density of states
 */
class DensityOfStates : public PhysicalModel
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


    //! Set the effective DOS
    void set_effective_DOS(double DOS);


    //! Get the total state density
    /*!
     * This method returns the total density of states, i.e.
     * the integral over all energies of the DOS.
     */
    double get_total_state_density(void);


    //! Get occupied states
    /*!
     * \return the density of occupied states in cm^-3
     * for the given quasi Fermi-level \c Ef and the derivative
     * w.r.t the quasi Fermi-level
     *
     * \param Ef quasi Fermi-level in eV
     * \param Epot potential energy (electrostatic energy) in eV
     * \param kT temperature in eV
     */
    std::pair<double, double>
    get_occupied_density_and_derivative(double Ef, double Epot, double kT,
        const Elem* elem, const Point& p, double kTlattice = -1) const;

    //! Get the occupied states
    /*!
     * In this version, the result is passed back in a vector. The size of the vector
     * will decide, what data is returned. The first element contains the density,
     * the second the first derivative, and the third the second derivative.
     * All derivatives are w.r.t the quasi Fermi-level
     */
    void get_occupied_density_and_derivative(std::vector<double>& result,
        double Ef, double Epot, double kT,
        const Elem* elem, const Point& p, double kTlattice = -1) const;


    //overloading for Trap.C
    //std::pair<double, double>
    //get_occupied_density_and_derivative(double Ef, double Epot,
    //    double kT, double kTlattice = -1) const;

    //! Get the particle name
    char get_particle(void) const;

    //! Get the particle spin
    double get_spin(void) const;


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

    //! Read the _fixed_DOS flag
    bool get_fixed_dos(void) const { return _fixed_DOS; }

    //! Get the thermoelectric power
    /*!
     * Must be called after \c calculate_density_and_derivative()
     */
    double get_thermoelectric_power(void) const;


  protected:

    //! Constructor
    DensityOfStates(const ModelOptions& options);

    //! A readable reference to the effective DOS parameter
    double& effective_dos(void) { return _effective_dos; }

    //! Get the effective DOS at given position
    double get_effective_dos(const Elem* elem, const Point& p) const;

    //! A RW reference to _fixed_DOS
    bool& fixed_dos(void) { return _fixed_DOS; }

    //! Calculate density and derivative
    virtual void
    calculate_density_and_derivative(std::vector<double>& result,
        double Ef, double Epot, double kT, double kTlattice,
        const Elem* elem, const Point& p) const = 0;

    //! overloading for Trap.C
    //virtual void
    //calculate_density_and_derivative(std::vector<double>& result,
    //    double Ef, double Epot, double kT, double kTlattice) const;

    //! Do we have a quantum density?
    bool& is_quantum_density(void);

    //! Set the thermoelectric power
    double& thermoelectric_power(void);

    //! Get or set the total state density
    /*!
     * DOS model \em must set the total state density!
     */
    double& total_state_density(void);

    //! The thermoelectric power
    /*!
     * The thermoelectric power is calculated in \c calculate_density_and_derivative()
     */
    mutable double _th_el_power;

  private:

    //! The reference energy
    std::vector<double> _reference_energy;

    //! The effective_mass
    std::vector<double> _effective_mass;

    //! The effective DOS
    double _effective_dos;

    //! A flag used to lock the effective DOS
    bool _fixed_DOS;

    //! We keep the name of the particle
    /*!
     * e -> electron
     * h -> hole
     */
    char _particle;

    //! The spin of the particle
    double _spin;


    //! Do or don't use quantum density
    bool _use_quantum;


    //! \c true if the last calculated density is a quantum density
    bool _is_quantum;


    //! The total density of states
    /*!
     * The density calculated from the module will always be limited
     * to this amount. In some cases this is artificial, however it is 
     * consistent with reality inasmuch as any band has a finite bandwidth
     * */
    double _total_density;

    //! A profile for the effective DOS
    ExternalProfile* _profile;

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

inline
void
DensityOfStates::set_effective_DOS(double DOS)
{
  if (!_fixed_DOS)
    _effective_dos = DOS;
}


inline
double
DensityOfStates::get_total_state_density(void)
{
  return _total_density;
}




inline
char
DensityOfStates::get_particle(void) const
{
  return _particle;
}

inline
double
DensityOfStates::get_spin(void) const
{
  return _spin;
}

inline
double&
DensityOfStates::total_state_density(void)
{
  return(_total_density);
}


inline
std::pair<double, double>
DensityOfStates::get_occupied_density_and_derivative(double Ef, double Epot,
    double kT, const Elem* elem, const Point& p, double kTlattice) const
{
  kTlattice = (kTlattice < 0) ? kT : kTlattice;
  std::vector<double> result(2, 0.0);
  get_occupied_density_and_derivative(result, Ef, Epot, kT, elem, p, kTlattice);

  return(std::make_pair(result[0], result[1]));
}

inline
void
DensityOfStates::get_occupied_density_and_derivative(std::vector<double>& result,
    double Ef, double Epot, double kT,
    const Elem* elem, const Point& p, double kTlattice) const
{
  kTlattice = (kTlattice < 0) ? kT : kTlattice;
  calculate_density_and_derivative(result, Ef, Epot, kT, kTlattice, elem, p);
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
  return(_is_quantum);
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


inline
double
DensityOfStates::get_thermoelectric_power(void) const
{
  return _th_el_power;
}

inline
double&
DensityOfStates::thermoelectric_power(void)
{
  return _th_el_power;
}


#endif // TC_DENSITYOFSTATES_H
