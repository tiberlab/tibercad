// $Id$

#ifndef _GENERICRECOMBINATION_H_
#define _GENERICRECOMBINATION_H_

#include "RecombinationModelInterface.h"
#include "TypeDefs.h"

class SimulationInterface;

//! Implementation of a generic recombination model with constant parameter
/*!
 * This class implements recombination processes that can be
 * modeled by \f[R=C\prod_i n_i^{\alpha_i}(1-1/kT\sum_i\pm\alpha_i\phi_i\f]
 */
class TBDLLOCAL GenericRecombination : public RecombinationModelInterface
{

  public:

    //! Destructor
    virtual ~GenericRecombination(void) {};

    //! Create a ConstantMobility object
    static GenericRecombination* create(const ModelOptions& options);

    
  protected:

    //! Constructor
    GenericRecombination(const ModelOptions& options);


    //! \copydoc RecombinationModelInterface::do_init()
    virtual void do_init(void) override;


    //! \copydoc RecombinationModelInterface::calculate_rate_and_derivatives()
    virtual void calculate_rate_and_derivatives(std::vector<double>& R,
        std::vector<std::vector<double>>& dPotentials) override;

  private:

    //! Recombination rate parameter
    double C_;

};



//
// inline methods
// 

inline
GenericRecombination::GenericRecombination(const ModelOptions& options)
  : RecombinationModelInterface(options),
    C_(0.0)
{
}


inline
GenericRecombination*
GenericRecombination::create(const ModelOptions& options)
{
  return new GenericRecombination(options);
}






#endif // _GENERICRECOMBINATION_H__
