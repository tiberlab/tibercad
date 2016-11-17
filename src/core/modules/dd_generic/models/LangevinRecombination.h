#ifndef _LANGEVINRECOMBINATION_H_
#define _LANGEVINRECOMBINATION_H_

#include "RecombinationModelInterface.h"
#include "TypeDefs.h"

class SimulationInterface;

//! Implementation of Langevin recombination
/*!
 * This class implements direct recombination processes that can be
 * modeled by \f[R_{Langevin}=\gamma \frac{2}{\varepsilon_r \varepsilon_0}(\mu_e + \mu_h)np\f]
 */
class TBDLLOCAL LangevinRecombination : public RecombinationModelInterface
{

  public:

    //! Destructor
    virtual ~LangevinRecombination(void) {};

    //! Create a LangevinRecombination object
    static LangevinRecombination* create(const ModelOptions& options);

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
    LangevinRecombination(const ModelOptions& options);

    //! \copydoc RecombinationModelInterface::read_database()
    virtual void read_database(void);

    //! \copydoc RecombinationModelInterface::do_init()
    virtual void do_init(void);

    //! \copydoc RecombinationModelInterface::do_reinit()
    virtual void do_reinit(void);

  private:

    typedef std::map<std::pair<SimulationInterface*, SimulationInterface*>,
        std::pair<unsigned int, double> > QRecMap;

    //! Relative permittivity from database
    double _er;

    //! gamma factor
    double _gamma;

};



//
// inline methods
// 

inline
LangevinRecombination::LangevinRecombination(const ModelOptions& options)
  : RecombinationModelInterface(options),
    _er(0.0),
    _gamma(1.0)
{
}


inline
LangevinRecombination*
LangevinRecombination::create(const ModelOptions& options)
{
  return new LangevinRecombination(options);
}

#endif // _LANGEVINRECOMBINATION_H__
