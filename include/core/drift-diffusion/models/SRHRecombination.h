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

    //! Destructor
    virtual ~SRHRecombination(void) {};

    //! Create a ConstantMobility object
    static SRHRecombination* create(void);

    //! \copydoc RecombinationModelInterface::get_net_recombination_rates()
    void get_net_recombination_rates(double& recomb_e, double& recomb_h);

    //! \copydoc
    //! RecombinationModelInterface::get_net_recombination_rate_derivatives()
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
    SRHRecombination(void);

    //! \copydoc RecombinationModelInterface::do_init()
    virtual void do_init(void);

    //! \copydoc RecombinationModelInterface::create_new()
    virtual PhysicalModelInterface* create_new(void) const;

    //! \copydoc RecombinationModelInterface::copy_from()
    virtual void copy_from(const PhysicalModelInterface* rhs);

    /*! \copydoc RecombinationModelInterface::calculate_VCA() */
    virtual void calculate_VCA(const PhysicalModelInterface* comp_A,
        const PhysicalModelInterface* comp_B, double xa);

    
  private:

    //! electron lifetime
    double _tau_n;

    //! hole lifetime
    double _tau_p;

};


//
// inline methods
// 


inline
SRHRecombination::SRHRecombination(void)
  : _tau_n(1e-9),
    _tau_p(1e-9)
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
  _tau_n = tau_n;
  _tau_p = tau_p;
}


inline
PhysicalModelInterface*
SRHRecombination::create_new(void) const
{
  return new SRHRecombination();
}


inline
void
SRHRecombination::copy_from(const PhysicalModelInterface* rhs)
{
  RecombinationModelInterface::copy_from(rhs);
  
  const SRHRecombination* mod = dynamic_cast<const SRHRecombination*>(rhs);
  _tau_n = mod->_tau_n;
  _tau_p = mod->_tau_p;
}


#endif // _SRHRECOMBINATION_H_
