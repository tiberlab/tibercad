#ifndef _HOSTGUESTRECOMBINATION_
#define _HOSTGUESTRECOMBINATION_

#include "RecombinationModelInterface.h"
#include "TypeDefs.h"

class SimulationInterface;

//! Implementation of HOST-GUEST recombination

class TBDLLOCAL HostGuestRecombination : public RecombinationModelInterface
{

  public:

    //! Destructor
    virtual ~HostGuestRecombination(void) {};

    //! Create a HostGuestRecombination object
    static HostGuestRecombination* create(const ModelOptions& options);

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
    HostGuestRecombination(const ModelOptions& options);

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

    //! sweep parameter for adiabatic switching
    double _alpha;

    //! Flag for adiabatic switching
    bool _adiabatic;

    //! Old solution for adiabatic switching;
    SimulationInterface* _saved_old_sim;

    //! Langevin factor;
    bool _langevin;


    //! The coupled simulation to use
    SimulationInterface* _coupled_sim;

    //! The current simulation
    SimulationInterface* _this_sim;

    //! Flag preventing infinite recursion
    static bool _coupled;

    //! Solution IDs
    ID _eDensity;
    ID _hDensity;
    ID _eDensity0;
    ID _hDensity0;
    ID _eMobility;
    ID _hMobility;
    ID _old_erec;
    ID _old_hrec;

};



//
// inline methods
// 

inline
HostGuestRecombination::HostGuestRecombination(const ModelOptions& options)
  : RecombinationModelInterface(options),
    _er(1.0),
    _gamma(1.0),
    _alpha(1.0),
    _langevin(true),
    _adiabatic(false)
{
}


inline
HostGuestRecombination*
HostGuestRecombination::create(const ModelOptions& options)
{
  return new HostGuestRecombination(options);
}

#endif // _HostGuestRecombination_H__
