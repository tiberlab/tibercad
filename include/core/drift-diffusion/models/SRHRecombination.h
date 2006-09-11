// $Id$

#ifndef _SRHRECOMBINATION_H_
#define _SRHRECOMBINATION_H_

#include "RecombinationModelInterface.h"
#include "TypeDefs.h"


//! Implementation of SRH recombination
/*!
 * This class implements Shockley-Read-Hall recombination processes that can be
 * modeled by 
 * \f[R_{SRH}=\frac{np - n_i^2}{(n+n_i)\tau_p + (p+n_i)\tau_n}\f]
 */
class SRHRecombination : public RecombinationModelInterface
{

  public:

    //! Constructor
    SRHRecombination(void);

    //! Destructor
    virtual ~SRHRecombination(void) {};

    //! \copydoc RecombinationModelInterface::get_net_recombination_rates()
    void get_net_recombination_rates(double& recomb_e, double& recomb_h);

    //! \copydoc
    //! RecombinationModelInterface::get_net_recombination_rate_derivatives()
    void get_net_recombination_rate_derivatives(
        std::vector<double>& recomb_e, std::vector<double>& recomb_h);

    //! Set the carrier lifetimes
    /*!
     * \param tau_n electron lifetime
     * \param tau_p hole lifetime
     */
    void set_SRH_parameters(double tau_n, double tau_p);

    //! \copydoc RecombinationModelInterface::set_model_options()
    virtual void set_model_options(const ModelOptions& options);

    //! \copydoc RecombinationModelInterface::get_name()
    virtual const std::string get_name(void) const;
    
  private:

    //! electron lifetime
    double _tau_n;

    //! hole lifetime
    double _tau_p;

};


inline
void
SRHRecombination::set_SRH_parameters(double tau_n, double tau_p)
{
  _tau_n = tau_n;
  _tau_p = tau_p;
}

#endif // _SRHRECOMBINATION_H_
