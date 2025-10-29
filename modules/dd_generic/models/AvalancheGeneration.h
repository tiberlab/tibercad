// $Id: AvalancheGeneration.h 4145 2015-10-02 11:53:20Z maufder $

#ifndef _AVALANCHEGENERATION_H_
#define _AVALANCHEGENERATION_H_

#include "RecombinationModelInterface.h"
#include "TypeDefs.h"


//! Implementation of Impact Ionization model
/*!
 * This class implements the impact ionization model processes according to
 *
 * \f{eqnarray*}
 * G_{II} & =& \sum_i\alpha_i|j_i| \\
 * \alpha_i & = & \gamma a_i e^{-\frac{\gamma b_i}{|E|}} \\
 * \f}
 *
 */
class TBDLLOCAL AvalancheGeneration : public RecombinationModelInterface
{

  public:

    //! Destructor
    virtual ~AvalancheGeneration(void);

    //! Create a ConstantMobility object
    static AvalancheGeneration* create(const ModelOptions& options);


  protected:

    //! Constructor
    AvalancheGeneration(const ModelOptions& options);

    //! \copydoc RecombinationModelInterface::read_database()
    virtual void read_database(void) override;

    //! \copydoc RecombinationModelInterface::do_init()
    virtual void do_init(void) override;


    //! \copydoc RecombinationModelInterface::calculate_rate_and_derivatives()
    virtual void calculate_rate_and_derivatives(std::vector<double>& R, std::vector<std::vector<double>>& dPotentials) override;


  private:

    //! Parameter \c a
    std::vector<double> _a_param;

    //! Parameter \c b
    std::vector<double> _b_param;

    //! Phonon energy
    double _w0;

};


//
// inline methods
//





inline
AvalancheGeneration*
AvalancheGeneration::create(const ModelOptions& options)
{
  return new AvalancheGeneration(options);
}




#endif // _AVALANCHEGENERATION_H_
