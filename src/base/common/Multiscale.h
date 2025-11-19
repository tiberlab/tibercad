// $Id$

#ifndef _MULTISCALE_H_
#define _MULTISCALE_H_


#include "tibercad/base/ModelOptions.h"
#include "tibercad/base/TypeDefs.h"

#include <set>

class SimulationInterface;

//! Class for the handling of multiscale couplings
class Multiscale
{

  public:

    //! The coupling method
    enum Method
    {
      UNKNOWN, //!< invalid method
      OVERLAP, //!< overlap method
      BRIDGE   //!< bridge method
    };

    //! Constructor
    explicit Multiscale(const ModelOptions& options);

    //! Reinit the multiscale coupling
    /*!
     * Can be called before solving \c sim.
     * If \c sim is one of the multiscale coupled models,
     * actions will be taken according to multiscale method.
     */
    void reinit(SimulationInterface* sim);


  private:

    //! The options
    ModelOptions _options;

    //! The coupling method
    Method _method;

    //! The macroscale model
    SimulationInterface* _macromodel;

    //! The microscale model
    SimulationInterface* _micromodel;

    //! The restriction of the macro domain
    std::set<ID> _restricted_macro_domain;

    //! The restricted variables
    std::vector<std::string> _restricted_variables;

    //! Initialize
    void init(void);

};


#endif // _MULTISCALE_H_
