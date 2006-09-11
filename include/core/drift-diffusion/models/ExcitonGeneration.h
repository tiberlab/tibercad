// $Id$

#ifndef _EXCITONGENERATION_H_
#define _EXCITONGENERATION_H_

#include "RecombinationModelInterface.h"
#include "TypeDefs.h"

class ExcitonTransport;


//! Implementation of Exciton generation
/*!
 * This class implements Exciton generation processes that can be
 * modeled by \f[G_{x}=Cnp\f]
 */
class ExcitonGeneration : public RecombinationModelInterface
{

  public:

    //! Constructor
    ExcitonGeneration(void);

    //! Destructor
    virtual ~ExcitonGeneration(void) {};

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
ExcitonGeneration::set_parameters(double C)
{
  _C = C;
}


#endif // _EXCITONGENERATION_H_
