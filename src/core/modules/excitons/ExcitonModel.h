// $Id$

#ifndef _EXCITONMODEL_H_
#define _EXCITONMODEL_H_

#include "SimulationOptions.h"
#include "ExcitonProperties.h"
#include "TiberCad.h"
#include "TypeDefs.h"
#include "ExcPolProps.h"

#include <string>
#include <map>

class SimulationInterface;

//! A simple Exciton model
class ExcitonModel : public ExcitonProperties
{
    
  public:
  
    //! The empty constructor.
    ExcitonModel(const ModelOptions& options);
       
    //! A default (empty) destructor.
    virtual ~ExcitonModel(void);

    //! This method creates an ExcitonModel object
    static ExcitonModel* create(const ModelOptions& options);

  protected:

    /*! \copydoc ExcitonProperties::prepare_element_data() */
    virtual void prepare_element_data(void);

    /*! \copydoc PhysicalModel::create_new() */
    virtual PhysicalModelInterface* create_new(void) const;

    /*! \copydoc ExcitonProperties::read_database() */
    virtual void read_database(void);

    /*! \copydoc ExcitonProperties::do_init() */
    virtual void do_init(void);

    virtual void do_recombination(void);

    virtual void do_mobility(void);

    
    //! Get the nonradiative recombination rate
    virtual double get_nonradiative_recombination_rate(void);

    //! Get the radiative recombination rate
    virtual double get_radiative_recombination_rate(void);
    
    //! Get the dissociation rate
    virtual double get_dissociation_rate(void);

    virtual double get_generation_rate();

    //! Get the dissociation rate
    virtual double get_exc_photon_scattering();

    //! Get the dissociation rate
    virtual double get_exc_exc_scattering();

    virtual double get_real_density() const;

    virtual double get_real_net_recombination_rate() const;

  private:
    ExcPolProps* excpolprops;


    //! Exciton radiative recombination time
    double _t_r;
    
    //! Exciton non-radiative recombination time
    double _t_nr;
    
    //! Exciton dissociation time
    double _t_diss;

    //! Exciton binding energy
    double _R;

    //! Exciton effective mass
    double _m;

    //! Exciton mobility
    double _mu;

    //! The ID for the generation model to be used
    ID _gen_model;

    //! The ID for the band gap variable
    ID _Eg_id;

    //! The DriftDiffusion simulation to be used
    SimulationInterface* _dd_sim;

    //! The band gap for each element
    std::map<const Elem*, double> _bandgap_data;


};


//
// inline members
//

inline
ExcitonModel::~ExcitonModel(void)
{
}


inline
PhysicalModelInterface*
ExcitonModel::create_new(void) const
{
  return new ExcitonModel(get_options());
}


inline
ExcitonModel*
ExcitonModel::create(const ModelOptions& options)
{
  return new ExcitonModel(options);
}


#endif /* _EXCITONMODEL_H_*/
