// $Id$

#ifndef _SRHRECOMBINATION_H_
#define _SRHRECOMBINATION_H_

#include "RecombinationModelInterface.h"
#include "TypeDefs.h"


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
 */
class TBDLEXPORT SRHRecombination : public RecombinationModelInterface
{

  public:

    //! Destructor
    virtual ~SRHRecombination(void) {};

    //! Create a ConstantMobility object
    static SRHRecombination* create(const ModelOptions& options);

    /*! \copydoc RecombinationModelInterface::get_net_recombination_rates() */
    void get_net_recombination_rates(double& recomb_e, double& recomb_h);

    /*! \copydoc
     * RecombinationModelInterface::get_net_recombination_rate_derivatives()
     */
    void get_net_recombination_rate_derivatives(
        std::vector<double>& recomb_e, std::vector<double>& recomb_h);

    //! Set the carrier lifetimes
    /*!
     * \deprecated { parameters will be read only from databeas or
     * ModelParameters }
     * \param tau_n electron lifetime
     * \param tau_p hole lifetime
     */
    void set_SRH_parameters(double tau_n, double tau_p);

    
  protected:

    //! Constructor
    TBDLLOCAL SRHRecombination(const ModelOptions& options);

    //! \copydoc RecombinationModelInterface::read_database()
    virtual void read_database(void);

    //! \copydoc RecombinationModelInterface::do_init()
    virtual void do_init(void);

    //! \copydoc RecombinationModelInterface::create_new()
    virtual PhysicalModelInterface* create_new(void) const;

    
  private:

    //! The electron lifetime
    double _tau_n;

    //! The hole lifetime
    double _tau_p;

    //! The trap level (from midgap)
    double _E_t;

    //! Temperature coefficient for temperature dependence, electrons
    double _Talpha_e;

    //! Temperature coefficient for temperature dependence, holes
    double _Talpha_h;

    //! Temperature coefficient for exponential temp. dependence, electrons
    double _Tcoeff_e;

    //! Temperature coefficient for exponential temp. dependence, holes
    double _Tcoeff_h;

};


//
// inline methods
// 


inline
SRHRecombination::SRHRecombination(const ModelOptions& options)
  : RecombinationModelInterface(options),
    _tau_n(1e-9),
    _tau_p(1e-9),
    _E_t(0.0),
    _Talpha_e(0.0),
    _Talpha_h(0.0),
    _Tcoeff_e(0.0),
    _Tcoeff_h(0.0)
{
}


inline
SRHRecombination*
SRHRecombination::create(const ModelOptions& options)
{
  return new SRHRecombination(options);
}


inline
void
SRHRecombination::set_SRH_parameters(double tau_n, double tau_p)
{
  _tau_n = tau_n;
  _tau_p = tau_p;
}


inline
PhysicalModelInterface*
SRHRecombination::create_new(void) const
{
  return new SRHRecombination(get_options());
}


#endif // _SRHRECOMBINATION_H_
