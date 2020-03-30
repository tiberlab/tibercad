// $Id: SRHRecombination.h 4145 2015-10-02 11:53:20Z maufder $

#ifndef _SRHRECOMBINATION_H_
#define _SRHRECOMBINATION_H_

#include "RecombinationModelInterface.h"
#include "TypeDefs.h"

class DensityOfStates;
class ExternalProfile;

//! Implementation of SRH recombination
/*!
 * This class implements Shockley-Read-Hall recombination processes that can be
 * modeled by
 * \f{eqnarray*}
 * R_{SRH} & =& \frac{np - n_i^2}{(n+n_i e^{(E_t-E_i)/k_BT})\tau_p +
 * (p+n_i e^{(E_i-E_t)/k_BT})\tau_n} \\
 * \tau_{n,p} & = & \tau_{n,p}^0 \left\{
 *     \begin{array}{l}
 *      \left(\frac{T}{T_0}\right)^{T_\alpha} \\
 *      e^{T_c(T/T_0 - 1)}
 *     \end{array}\right.
 * \f}
 * The recombination times are doping dependent, e.g.:
 * \f[
 * \tau_n^0 = \tau_{min,n} + \frac{\tau_{max,n} -
 *         \tau_{min,n}}{1 + (N/N_{ref})^\gamma}
 * \f]
 *
 * Trap-assisted tunneling is implemented following Hurkx et al., "A New
 *   Recombination Model for Device Simulation Including Tunneling",
 *   IEEE Trans. on Electron Devices, 39, 331-338, 1992
 */
class TBDLLOCAL SRHRecombination : public RecombinationModelInterface
{

  public:

    //! Destructor
    virtual ~SRHRecombination(void);

    //! Create a ConstantMobility object
    static SRHRecombination* create(const ModelOptions& options);


  protected:

    //! Constructor
    SRHRecombination(const ModelOptions& options);

    //! \copydoc RecombinationModelInterface::read_database()
    virtual void read_database(void) override;

    //! \copydoc RecombinationModelInterface::read_interface_database()
    virtual void read_interface_database(void) override;

    //! \copydoc RecombinationModelInterface::prepare_submodels()
    virtual void prepare_submodels(void) override;

    //! \copydoc RecombinationModelInterface::do_init()
    virtual void do_init(void) override;

    //! \copydoc RecombinationModelInterface::do_init_interface()
    virtual void do_init_interface(const Material* comp_A,
            const Material* comp_B) override;

    //! We do the VCA on the recombination times and not on the parameters
    virtual void do_init_alloy(const PhysicalModel* comp_A,
        const PhysicalModel* comp_B, double xa) override;

    //! \copydoc RecombinationModelInterface::calculate_rate_and_derivatives()
    virtual void calculate_rate_and_derivatives(std::vector<double>& R, std::vector<std::vector<double>>& dPotentials) override;


  private:

    struct TrapAssisted
    {
      TrapAssisted(void);

      //! TAT tunneling mass
      double m_trap;

      //! calculate Gamma
      /*!
       * F is electric field in V/m, T temperature in K,
       * Et trap level from conduction band
       */
      double get_gamma(double F, double T, double Et);
    };

    //! True if this is a trap generated model
    bool _trap;

    //! The electron lifetime
    double _tau_n;

    //! The hole lifetime
    double _tau_p;

    //! electron capture cross section
    double _sigma_n;

    //! hole capture cross section
    double _sigma_p;

    //! The trap level
    double _E_t;

    //! The trap density
    double _density;

    //! Trap-to-CB generation
    double _gen_TC;

    //! VB-to-trap generation
    double _gen_VT;


    //! The energy reference
    char _energy_reference;

    //! Temperature coefficient for temperature dependence, electrons
    double _Talpha_e;

    //! Temperature coefficient for temperature dependence, holes
    double _Talpha_h;

    //! Temperature coefficient for exponential temp. dependence, electrons
    double _Tcoeff_e;

    //! Temperature coefficient for exponential temp. dependence, holes
    double _Tcoeff_h;

    //! The TAT coefficients, if present
    TrapAssisted* _tat;

    //! The density of states, if provided
    DensityOfStates* _dos;

    //! Trap density profile
    ExternalProfile* _profile;

    //! Get the trap level
    double get_trap_level(void);
};


//
// inline methods
//





inline
SRHRecombination*
SRHRecombination::create(const ModelOptions& options)
{
  return new SRHRecombination(options);
}




#endif // _SRHRECOMBINATION_H_
