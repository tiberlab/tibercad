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

    //! Create a ConstantMobility object
    static ExcitonGeneration* create(void);

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

    //! \copydoc RecombinationModelInterface::do_init()
    virtual void do_init(void);

    //! \copydoc RecombinationModelInterface::create_new()
    virtual PhysicalModelInterface* create_new(void) const;

    /*! \copydoc RecombinationModelInterface::do_init_alloy() */
    virtual void do_init_alloy(const PhysicalModelInterface* comp_A,
        const PhysicalModelInterface* comp_B, double xa);

    
  private:

    //! Recombination rate parameter
    double C_;

};


//
// inline methods
// 


inline
ExcitonGeneration::ExcitonGeneration(void)
  : C_(1e-10)
{
}


inline
ExcitonGeneration*
ExcitonGeneration::create(void)
{
  return new ExcitonGeneration();
}


inline
void
ExcitonGeneration::set_parameters(double C)
{
  C_ = C;
}


inline
PhysicalModelInterface*
ExcitonGeneration::create_new(void) const
{
  return new ExcitonGeneration();
}


#endif // _EXCITONGENERATION_H_
