// $Id$

#ifndef _OPTICALGENERATION_H_
#define _OPTICALGENERATION_H_

#include "RecombinationModelInterface.h"
#include "Variable.h"
#include "TypeDefs.h"


//! Implementation of optical generation
/*!
 * This class implements optical generation processes that can be
 * modeled by \f[G_{x}= G]
 *
 * It is derived from Variable to be able to make a sweep over the 
 * generation rate.
 */
class OpticalGeneration : public RecombinationModelInterface, public Variable
{

  public:

    //! Constructor
    OpticalGeneration(void);

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


  protected:

    //! \copydoc RecombinationModelInterface::do_init()
    virtual void do_init(void);

    //! \copydoc RecombinationModelInterface::create_new()
    virtual PhysicalModelInterface* create_new(void) const;

    /*! \copydoc RecombinationModelInterface::do_init_alloy() */
    virtual void do_init_alloy(const PhysicalModelInterface* comp_A,
        const PhysicalModelInterface* comp_B, double xa);


    /*! \copydoc Variable::set_variable_value() */
    virtual void set_variable_value(double value, ID id = 0);


    /*! \copydoc Variable::set_variable_value() */
    virtual double get_variable_value(ID id = 0);


  private:

    //! Generation rate parameter
    double G_;

};



//
// inline methods
// 

inline
OpticalGeneration::OpticalGeneration(void)
  : G_(0.0)
{
}


inline
OpticalGeneration*
OpticalGeneration::create(void)
{
  return new OpticalGeneration();
}



inline
PhysicalModelInterface*
OpticalGeneration::create_new(void) const
{
  return new OpticalGeneration();
}


inline
void
OpticalGeneration::set_variable_value(double value, ID id)
{
  ignore_unused_variable(id);
  G_ = value;
}


inline
double
OpticalGeneration::get_variable_value(ID id)
{
  ignore_unused_variable(id);
  return G_;
}


#endif // _OPTICALGENERATION_H_
