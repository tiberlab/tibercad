// $Id$

#ifndef _OPTICALGENERATION_H_
#define _OPTICALGENERATION_H_

#include "RecombinationModelInterface.h"
#include "TypeDefs.h"
#include "vector_value.h"
#include "SimulationInterface.h"
#include "point.h"

// forward declarations
class Elem;



//! Implementation of optical generation
/*!
 * This class implements optical generation processes that can be
 * modeled by \f[G_{x}= G]
 *
 */
class TBDLLOCAL OpticalGeneration : public RecombinationModelInterface
{

  public:

    //! Destructor
    virtual ~OpticalGeneration(void) {};

    //! Create a ConstantMobility object
    static OpticalGeneration* create(const ModelOptions& options);

    //! \copydoc RecombinationModelInterface::get_net_recombination_rates()
    void get_net_recombination_rates(double& recomb_e, double& recomb_h);

    /*!
     * \copydoc
     * RecombinationModelInterface::get_net_recombination_rate_derivatives()
     */
    void get_net_recombination_rate_derivatives(
        std::vector<double>& recomb_e, std::vector<double>& recomb_h);

    //! File to read the external generation
//    double read_file(char *filename, double p);
     double read_file(void);

  protected:

    //! Constructor
    OpticalGeneration(const ModelOptions& options);

    //! \copydoc RecombinationModelInterface::do_init()
    virtual void do_init(void);


  private:

    //! Generation rate parameter
    double _generation;

    //! Flag to decide to upload an external generation file
    bool _read_file;

    //! Number of suns
    double _sun;

    //! The generation model
    std::vector<SimulationInterface*> _generation_model;

    //! The solution ID of the generation models variable
    std::vector<ID> _gen_id;


};



//
// inline methods
//

inline
OpticalGeneration::OpticalGeneration(const ModelOptions& options)
  : RecombinationModelInterface(options),
    _generation(0.0),
    _sun(0.0),
    _read_file(false)
{
}


inline
OpticalGeneration*
OpticalGeneration::create(const ModelOptions& options)
{
  return new OpticalGeneration(options);
}





#endif // _OPTICALGENERATION_H_
