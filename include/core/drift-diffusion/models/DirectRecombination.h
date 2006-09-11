// $Id$

#ifndef _DIRECTRECOMBINATION_H_
#define _DIRECTRECOMBINATION_H_

#include "RecombinationModelInterface.h"
#include "TypeDefs.h"


//! Implementation of direct recombination
/*!
 * This class implements direct recombination processes that can be
 * modeled by \f[R_{direct}=C(np-n_i^2)\f]
 */
class DirectRecombination : public RecombinationModelInterface
{

  public:

    //! Constructor
    DirectRecombination(void);

    //! Destructor
    virtual ~DirectRecombination(void) {};

    //! \copydoc RecombinationModelInterface::get_net_recombination_rates()
    void get_net_recombination_rates(double& recomb_e, double& recomb_h);

    /*!
     * \copydoc
     * RecombinationModelInterface::get_net_recombination_rate_derivatives()
     */
    void get_net_recombination_rate_derivatives(
        std::vector<double>& recomb_e, std::vector<double>& recomb_h);

    //! Set the direct recombination parameters
    void set_parameters(double C);

    //! \copydoc RecombinationModelInterface::set_model_options()
    virtual void set_model_options(const ModelOptions& options);

    //! \copydoc RecombinationModelInterface::get_name()
    virtual const std::string get_name(void) const;
    
  private:

    //! Recombination rate parameter
    double _C;

};


inline
void
DirectRecombination::set_parameters(double C)
{
  _C = C;
}

#endif // _DIRECTRECOMBINATION_H_
