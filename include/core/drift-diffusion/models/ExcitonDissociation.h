// $Id$

#ifndef _EXCITONDISSOCIATION_H_
#define _EXCITONDISSOCIATION_H_

#include "RecombinationModelInterface.h"
#include "TypeDefs.h"

class ExcitonTransport;


//! Implementation of Exciton dissociation
/*!
 * This class implements Exciton dissociation process
 */
class ExcitonDissociation : public RecombinationModelInterface
{

  public:

    //! Constructor
    ExcitonDissociation(void);

    //! Destructor
    virtual ~ExcitonDissociation(void) {};

    //! \copydoc RecombinationModelInterface::get_net_recombination_rates()
    void get_net_recombination_rates(double& recomb_e, double& recomb_h);

    /*!
     * \copydoc
     * RecombinationModelInterface::get_net_recombination_rate_derivatives()
     */
    void get_net_recombination_rate_derivatives(
        std::vector<double>& recomb_e, std::vector<double>& recomb_h);

    //! Set the \c ExcitonTransport to be used
    void set_exciton_transport(ExcitonTransport* simulator);

    //! \copydoc RecombinationModelInterface::set_model_options()
    virtual void set_model_options(const ModelOptions& options);

    //! \copydoc RecombinationModelInterface::get_name()
    virtual const std::string get_name(void) const;

  private:

    //! Damping factor
    double _d;

    //! Trapping probability
    /*!
     * \f$\alpha\f$ is the percentage of dissociating excitons that
     * gets trapped instead of creating free electron hole pairs
     */
    double _a;

    //! The \c ExcitonTransport to use
    ExcitonTransport* _exciton_sim;

};



inline
void
ExcitonDissociation::set_exciton_transport(ExcitonTransport* simulator)
{
  _exciton_sim = simulator;
}


#endif // _EXCITONDISSOCIATION_H_
