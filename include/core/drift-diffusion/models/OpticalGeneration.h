// $Id: ExcitonGeneration.h 70 2006-07-14 16:57:35Z maufder $

#ifndef _OPTICALGENERATION_H_
#define _OPTICALENERATION_H_

#include "RecombinationModelInterface.h"
#include "TypeDefs.h"


//! Implementation of optical generation
/*!
 * This class implements optical generation processes that can be
 * modeled by \f[G_{x}= G]
 */
class OpticalGeneration : public RecombinationModelInterface
{

  public:

    //! Constructor
    OpticalGeneration(void);

    //! Destructor
    virtual ~OpticalGeneration(void) {};

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

    //! Generation rate parameter
    double _G;

};


inline
void
OpticalGeneration::set_parameters(double G)
{
  _G = G;
}

#endif // _OPTICALGENERATION_H_
