#ifndef _EXCITONTTAINTERPOLATED_H_
#define _EXCITONTTAINTERPOLATED_H_

#include "RecombinationModelInterface.h"
#include "TypeDefs.h"

class SimulationInterface;

//! Implementation of direct recombination
/*!
 * This class implements direct recombination processes that can be
 * modeled by \f[R_{direct}=C(np-n_i^2)\f]
 */
class TBDLLOCAL ExcitonTTAInterpolated : public RecombinationModelInterface
{

  public:

    //! Destructor
    virtual ~ExcitonTTAInterpolated(void) {};

    //! Create a ConstantMobility object
    static ExcitonTTAInterpolated* create(const ModelOptions& options);

    
  protected:

    //! Constructor
    ExcitonTTAInterpolated(const ModelOptions& options);

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

    //! The interpolation simulation to use
    SimulationInterface* _interpolation_sim;

    //! Model name
    std::string _model_name;

    //Variables
    std::string _temperature_var;
    std::string _density_var;

    ID _model_id;
    ID _temperature_id;
    ID _density_id;



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
ExcitonTTAInterpolated::ExcitonTTAInterpolated(const ModelOptions& options)
  : RecombinationModelInterface(options),
    _quantum_optics(NULL)
{
}


inline
ExcitonTTAInterpolated*
ExcitonTTAInterpolated::create(const ModelOptions& options)
{
  return new ExcitonTTAInterpolated(options);
}






#endif // _EXCITONTTAINTERPOLATED_H__
