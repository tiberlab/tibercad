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
class TBDLLOCAL ExcitonGeneration : public RecombinationModelInterface
{

  public:

    //! Destructor
    virtual ~ExcitonGeneration(void) {};

    //! Create a ConstantMobility object
    static ExcitonGeneration* create(const ModelOptions& options);

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


  protected:

    //! Constructor
    ExcitonGeneration(const ModelOptions& options);

    //! \copydoc RecombinationModelInterface::do_init()
    virtual void do_init(void);

    
  private:

    //! Recombination rate parameter
    double C_;

};


//
// inline methods
// 


inline
ExcitonGeneration::ExcitonGeneration(const ModelOptions& options)
  : RecombinationModelInterface(options),
    C_(1e-10)
{
}


inline
ExcitonGeneration*
ExcitonGeneration::create(const ModelOptions& options)
{
  return new ExcitonGeneration(options);
}


inline
void
ExcitonGeneration::set_parameters(double C)
{
  C_ = C;
}



#endif // _EXCITONGENERATION_H_
