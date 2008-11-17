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

    //! Constructor
    TBDLLOCAL SRHRecombination(void);

    //! Destructor
    virtual ~SRHRecombination(void) {};

    //! Create a ConstantMobility object
    static SRHRecombination* create(void);

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

    //! \copydoc RecombinationModelInterface::read_database()
    virtual void read_database(void);

    //! \copydoc RecombinationModelInterface::do_init()
    virtual void do_init(void);

    //! \copydoc RecombinationModelInterface::create_new()
    virtual PhysicalModelInterface* create_new(void) const;

    /*! \copydoc RecombinationModelInterface::do_init_alloy() */
    virtual void do_init_alloy(const PhysicalModelInterface* comp_A,
        const PhysicalModelInterface* comp_B, double xa);

    
  private:

    //! The electron lifetime
    double tau_n_;

    //! The hole lifetime
    double tau_p_;

    //! The trap level (from midgap)
    double E_t_;

    //! Temperature coefficient for temperature dependence, electrons
    double Talpha_e_;

    //! Temperature coefficient for temperature dependence, holes
    double Talpha_h_;

    //! Temperature coefficient for exponential temperature dependence
    //double Tcoeff_;

};


//
// inline methods
// 


inline
SRHRecombination::SRHRecombination(void)
  : tau_n_(1e-9),
    tau_p_(1e-9),
    E_t_(0.0),
    Talpha_e_(0.0),
    Talpha_h_(0.0)
    //_Tcoeeff(0.0)
{
}


inline
SRHRecombination*
SRHRecombination::create(void)
{
  return new SRHRecombination();
}


inline
void
SRHRecombination::set_SRH_parameters(double tau_n, double tau_p)
{
  tau_n_ = tau_n;
  tau_p_ = tau_p;
}


inline
PhysicalModelInterface*
SRHRecombination::create_new(void) const
{
  return new SRHRecombination();
}


#endif // _SRHRECOMBINATION_H_
