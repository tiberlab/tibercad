// $Id$

#ifndef _DIRECTRECOMBINATION_H_
#define _DIRECTRECOMBINATION_H_

#include "RecombinationModelInterface.h"
#include "TypeDefs.h"


//! Implementation of direct recombination
/*!
 * This class implements direct recombination processes that can be
 * modeled by \f[R_{direct}=C(np-n_i^2)\f]
 */
class DirectRecombination : public RecombinationModelInterface
{

  public:

    //! Constructor
    DirectRecombination(void);

    //! Destructor
    virtual ~DirectRecombination(void) {};

    //! Create a ConstantMobility object
    static DirectRecombination* create(void);

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

    //! \copydoc RecombinationModelInterface::read_database()
    virtual void read_database(void);

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

    //! Recombination rate parameter
    double C_;

};



//
// inline methods
// 

inline
DirectRecombination::DirectRecombination(void)
  : C_(0.0)
{
}


inline
DirectRecombination*
DirectRecombination::create(void)
{
  return new DirectRecombination();
}


inline
void
DirectRecombination::set_parameters(double C)
{
  C_ = C;
}


inline
PhysicalModelInterface*
DirectRecombination::create_new(void) const
{
  return new DirectRecombination();
}


inline
void
DirectRecombination::copy_from(const PhysicalModelInterface* rhs)
{
  RecombinationModelInterface::copy_from(rhs);
  
  const DirectRecombination* mod = dynamic_cast<const DirectRecombination*>(rhs);
  C_ = mod->C_;
}




#endif // _DIRECTRECOMBINATION_H__
