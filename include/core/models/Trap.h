// $Id$

#ifndef _TRAP_H_
#define _TRAP_H_


#include "PhysicalModel.h"
#include "SimulationInterface.h"
#include "ExternalProfile.h"

class DensityOfStates;
class Particle;

/*!
 * \brief Base class for traps in semiconductors
 *
 * The implemented traps in this model are:
 * \li neutral electron traps, type = ``eNeutral''
 * \li neutral hole traps, type = ``hNeutral''
 * \li donor like traps, type = ``donor''
 * \li acceptor like traps, type = ``acceptor''
 *
 * The trap level is referred to midgap as default.
 */
class TBDLEXPORT Trap : public PhysicalModel
{

  public:

    //! The trap type
    enum Type
    {
      NEUTRAL,  //!< neutral if unoccupied
      CHARGED,  //!< charged if unoccupied
      FIXED     //!< a fixed charge density
    };

    //! Destructor
    virtual ~Trap(void);


    //! Creator function
    static Trap* create(const ModelOptions& options);


    /*!
     * \brief Set the relevant band edge energies
     *
     * \param Ec the conduction band edge in eV
     * \param Ev the valence band edge in eV
     * \param phi the electrostatic potential
     */
    void set_energies(double Ec, double Ev, double phi);


    //! Get ionized density (= charge density, cm^-3)
    /*!
     * \return the charge density in cm^-3 including sign
     * \param elem the mesh element
     * \param p the point inside \c elem
     * \param el the electron population
     * \param hl the hole population
     * \param will be filled with the partial derivatives: [d/dn d/dp d/dEfn d/dEfp]
     */
    double get_ionized_density_and_derivative(const Elem* elem, const Point& p,
        const Particle& el, const Particle& hl,
        std::vector<double>& derivatives) const;

    //! Get the particle type
    char get_particle(void) const;


    //! Get the trap type
    Type get_type(void) const;


  protected:

    //! Constructor
    Trap(const ModelOptions& options);


    //! Initialize
    virtual void do_init(void);

    //! Create DOS model (if present)
    virtual void prepare_submodels(void);


    //! Set the trap type
    void set_type(Type type);


  private:

    //! The trap density in cm^-3
    double _density;

    //! The trap profile function
    ExternalProfile* _profile;

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

    //! The electrostatic potential
    double _phi;


    //! e cross section
    double _sigma_n;

    //! h cross section
    double _sigma_p;

    //! e capture rate
    double _tau_n;

    //! e capture rate
    double _tau_p;

    //! e thermal velocity
    double _e_vth;

    //! h thermal velocity
    double _h_vth;


    //! Trap-to-CB generation
    double _gen_TC;

    //! VB-to-trap generation
    double _gen_VT;


    //! The density of states
    DensityOfStates* _dos;


    //! Calculate the trap level
    double _trap_level(void) const TBDLLOCAL;

    //! The simulation providing external electron and hole densities
    SimulationInterface* _ext_dens_sim;

    //! Flag preventing infinite recursion
    static bool _coupled;

    //! Solution IDs
    ID _eDensity;
    ID _hDensity;

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
Trap::set_energies(double Ec, double Ev, double phi)
{
  _Ec = Ec;
  _Ev = Ev;
  _phi = phi;

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
