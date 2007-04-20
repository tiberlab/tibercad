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

    //! Create a ConstantMobility object
    static ExcitonDissociation* create(void);

    //! \copydoc RecombinationModelInterface::get_net_recombination_rates()
    virtual void get_net_recombination_rates(double& recomb_e, double& recomb_h);

    /*!
     * \copydoc
     * RecombinationModelInterface::get_net_recombination_rate_derivatives()
     */
    virtual void get_net_recombination_rate_derivatives(
        std::vector<double>& recomb_e, std::vector<double>& recomb_h);

    
  protected:
    
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

    //! Damping factor
    double d_;

    //! The \c ExcitonTransport to use
    ExcitonTransport* exciton_sim_;

};



//
// inline methods
// 

inline
ExcitonDissociation::ExcitonDissociation(void)
  : d_(1.0),
    exciton_sim_(NULL)
{
}

inline
ExcitonDissociation*
ExcitonDissociation::create(void)
{
  return new ExcitonDissociation();
}



inline
PhysicalModelInterface*
ExcitonDissociation::create_new(void) const
{
  return new ExcitonDissociation();
}


inline
void
ExcitonDissociation::copy_from(const PhysicalModelInterface* rhs)
{
  RecombinationModelInterface::copy_from(rhs);
  
  const ExcitonDissociation* mod = dynamic_cast<const ExcitonDissociation*>(rhs);
  exciton_sim_ = mod->exciton_sim_;
  d_ = mod->d_;
}


#endif // _EXCITONDISSOCIATION_H_
