// $Id$

#ifndef _OPTICALGENERATION_H_
#define _OPTICALGENERATION_H_

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

    //! Destructor
    virtual ~OpticalGeneration(void) {};

    //! Create a ConstantMobility object
    static OpticalGeneration* create(void);

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
    OpticalGeneration(void);

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

    //! Generation rate parameter
    double _G;

};



//
// inline methods
// 

inline
OpticalGeneration::OpticalGeneration(void)
  : _G(1e-10)
{
}


inline
OpticalGeneration*
OpticalGeneration::create(void)
{
  return new OpticalGeneration();
}


inline
void
OpticalGeneration::set_parameters(double G)
{
  _G = G;
}


inline
PhysicalModelInterface*
OpticalGeneration::create_new(void) const
{
  return new OpticalGeneration();
}


inline
void
OpticalGeneration::copy_from(const PhysicalModelInterface* rhs)
{
  RecombinationModelInterface::copy_from(rhs);
  
  const OpticalGeneration* mod = dynamic_cast<const OpticalGeneration*>(rhs);
  _G = mod->_G;
}


#endif // _OPTICALGENERATION_H_
