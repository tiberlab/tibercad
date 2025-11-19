// $Id: ExcitonGeneration.h 3414 2012-09-10 20:40:28Z maufder $

#ifndef _EXCITONGENERATION_H_
#define _EXCITONGENERATION_H_

#include "RecombinationModelInterface.h"
#include "tibercad/base/TypeDefs.h"

class SimulationInterface;

//! Implementation of exciton generation/dissociation from free carrier gas
/*!
 * This class implements generation and dissociation of excitons
 * modeled by \f[R=\gamma np(1 - exp(\frac{E_{Fx} - E_{Fn} + E_{Fp}}{k_BT}))\f]
 *
 * In the input file, the carriers have to be provided in the order
 * \c electron, \c hole, \c exciton
 */
class TBDLLOCAL ExcitonGeneration : public RecombinationModelInterface
{

  public:

    //! Destructor
    virtual ~ExcitonGeneration(void) {};

    //! Create a ConstantMobility object
    static ExcitonGeneration* create(const ModelOptions& options);

    
  protected:

    //! Constructor
    ExcitonGeneration(const ModelOptions& options);

    //! \copydoc RecombinationModelInterface::read_database()
    virtual void read_database(void) override;

    //! \copydoc RecombinationModelInterface::do_init()
    virtual void do_init(void) override;

    //! \copydoc RecombinationModelInterface::calculate_rate_and_derivatives()
    virtual void calculate_rate_and_derivatives(std::vector<double>& R,
        std::vector<std::vector<double>>& dPotentials) override;

  private:

    //! Recombination rate parameters
    double  _gamma;

    bool _stat_fac;

};



//
// inline methods
// 

inline
ExcitonGeneration::ExcitonGeneration(const ModelOptions& options)
  : RecombinationModelInterface(options),
    _stat_fac(true)
{
}


inline
ExcitonGeneration*
ExcitonGeneration::create(const ModelOptions& options)
{
  return new ExcitonGeneration(options);
}






#endif // _EXCITONGENERATION_H__
