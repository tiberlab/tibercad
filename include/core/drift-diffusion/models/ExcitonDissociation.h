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

    //! Destructor
    virtual ~ExcitonDissociation(void) {};

    //! Create a ConstantMobility object
    static ExcitonDissociation* create(void);

    //! \copydoc RecombinationModelInterface::get_net_recombination_rates()
    void get_net_recombination_rates(double& recomb_e, double& recomb_h);

    /*!
     * \copydoc
     * RecombinationModelInterface::get_net_recombination_rate_derivatives()
     */
    void get_net_recombination_rate_derivatives(
        std::vector<double>& recomb_e, std::vector<double>& recomb_h);

    
  protected:

    //! Constructor
    ExcitonDissociation(void);
    
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
    double _d;

    //! The \c ExcitonTransport to use
    ExcitonTransport* _exciton_sim;

};



//
// inline methods
// 

inline
ExcitonDissociation::ExcitonDissociation(void)
  : _d(1.0),
    _exciton_sim(NULL)
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
  _exciton_sim = mod->_exciton_sim;
  _d = mod->_d;
}


#endif // _EXCITONDISSOCIATION_H_
