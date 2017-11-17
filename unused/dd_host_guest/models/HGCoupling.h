// $Id: HGCoupling.h 3414 2012-09-10 20:40:28Z maufder $

#ifndef _HGCOUPLING_H_
#define _HGCOUPLING_H_

#include "RecombinationModelInterface.h"
#include "TypeDefs.h"


class TBDLLOCAL HGCoupling : public RecombinationModelInterface
{

  public:

    //! Destructor
    virtual ~HGCoupling(void) {};

    //! Create a ConstantMobility object
    static HGCoupling* create(const ModelOptions& options);

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
    HGCoupling(const ModelOptions& options);

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
HGCoupling::HGCoupling(const ModelOptions& options)
  : RecombinationModelInterface(options),
    C_(0.0)
{
  if (get_option("bands", "e") == std::string("h"))
    this->set_coupling("hHG");
  else
    this->set_coupling("eHG");
}


inline
HGCoupling*
HGCoupling::create(const ModelOptions& options)
{
  return new HGCoupling(options);
}






#endif // _DIRECTRECOMBINATION_H__
