// $Id: ExcitonISC.h 3414 2012-09-10 20:40:28Z maufder $

#ifndef _EXCITONISC_H_
#define _EXCITONISC_H_

#include "RecombinationModelInterface.h"
#include "TypeDefs.h"

class SimulationInterface;

//! Implementation of direct recombination
/*!
 * This class implements direct recombination processes that can be
 * modeled by \f[R_{direct}=C(np-n_i^2)\f]
 */
class TBDLLOCAL ExcitonISC : public RecombinationModelInterface
{

  public:

    //! Destructor
    virtual ~ExcitonISC(void) {};

    //! Create a ConstantMobility object
    static ExcitonISC* create(const ModelOptions& options);


  protected:

    //! Constructor
    ExcitonISC(const ModelOptions& options);

    //! \copydoc RecombinationModelInterface::read_database()
    virtual void read_database(void) override;

    //! \copydoc RecombinationModelInterface::do_init()
    virtual void do_init(void) override;

    //! \copydoc RecombinationModelInterface::do_reinit()
    virtual void do_reinit(void) override;

    //! \copydoc RecombinationModelInterface::calculate_rate_and_derivatives()
    virtual void calculate_rate_and_derivatives(std::vector<double>& R, std::vector<std::vector<double>>& dPotentials) override;

  private:

    typedef std::map<std::pair<SimulationInterface*, SimulationInterface*>,
        std::pair<unsigned int, double> > QRecMap;

    //! Recombination rate parameter
    double  _C;

    //! The quantum optics simulation, if available
    SimulationInterface* _quantum_optics;

    //! The solution ID for the optical recombination
    ID _rec_id;

    //! A static map to put quantum recombination in
    /*!
     * This map is used so as to not calculate the same quantity
     * several times.
     */
    static QRecMap _qrec_vals;

};



//
// inline methods
// 

inline
ExcitonISC::ExcitonISC(const ModelOptions& options)
  : RecombinationModelInterface(options),
    _quantum_optics(NULL)
{
}


inline
ExcitonISC*
ExcitonISC::create(const ModelOptions& options)
{
  return new ExcitonISC(options);
}






#endif // _EXCITONISC_H__
