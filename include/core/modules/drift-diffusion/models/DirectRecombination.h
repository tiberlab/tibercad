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
class TBDLLOCAL DirectRecombination : public RecombinationModelInterface
{

  public:

    //! Destructor
    virtual ~DirectRecombination(void) {};

    //! Create a ConstantMobility object
    static DirectRecombination* create(const ModelOptions& options);

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
    DirectRecombination(const ModelOptions& options);

    //! \copydoc RecombinationModelInterface::read_database()
    virtual void read_database(void);

    //! \copydoc RecombinationModelInterface::do_init()
    virtual void do_init(void);

    //! \copydoc RecombinationModelInterface::create_new()
    virtual PhysicalModelInterface* create_new(void) const;


  private:

    //! Recombination rate parameter
    double C_;

};



//
// inline methods
// 

inline
DirectRecombination::DirectRecombination(const ModelOptions& options)
  : RecombinationModelInterface(options),
    C_(0.0)
{
}


inline
DirectRecombination*
DirectRecombination::create(const ModelOptions& options)
{
  return new DirectRecombination(options);
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
  return new DirectRecombination(get_options());
}




#endif // _DIRECTRECOMBINATION_H__
