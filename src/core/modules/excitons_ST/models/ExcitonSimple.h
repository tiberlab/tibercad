// $Id: ExcitonSimple.h 4192 2015-12-10 11:11:18Z maufder $

#ifndef _EXCITONSIMPLE_H_
#define _EXCITONSIMPLE_H_

#include "SimulationOptions.h"
#include "ExcitonProperties.h"
#include "TiberCad.h"
#include "TypeDefs.h"

#include <string>
#include <map>

class SimulationInterface;

//! A simple Exciton model
class ExcitonSimple : public ExcitonProperties
{
    
  public:
  
    //! The empty constructor.
    ExcitonSimple(const ModelOptions& options);
       
    //! A default (empty) destructor.
    virtual ~ExcitonSimple(void);

    //! This method creates an ExcitonSimple object
    static ExcitonSimple* create(const ModelOptions& options);

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

    virtual void do_diffusion(void);

    
    //! Get the singlet nonradiative recombination rate
    virtual double get_s_nonradiative_recombination_rate(void);

    //! Get the singlet radiative recombination rate
    virtual double get_s_radiative_recombination_rate(void);
    
    //! Get the singlet dissociation rate
    virtual double get_s_dissociation_rate(void);

    //! Get the ISC rate
    virtual double get_isc_rate(void);

    //! Get the ISC rate derivative
    virtual double get_isc_rate_derivative(void);

    //! Get the singlet generation rate

    virtual double get_s_generation_rate();

    //! Get the triplet nonradiative recombination rate
    virtual double get_t_nonradiative_recombination_rate(void);

    //! Get the triplet radiative recombination rate
    virtual double get_t_radiative_recombination_rate(void);
    
    //! Get the triplet dissociation rate
    virtual double get_t_dissociation_rate(void);

    //! Get the triplet generation rate

    virtual double get_t_generation_rate();


  private:

    //! Singlet radiative recombination time
    double _ts_r;
    
    //! Singlet non-radiative recombination time
    double _ts_nr;
    
    //! Singlet dissociation time
    double _ts_diss;

    //! Triplet radiative recombination time
    double _tt_r;
    
    //! Triplet non-radiative recombination time
    double _tt_nr;
    
    //! Triplet dissociation time
    double _tt_diss;

    //! Exciton binding energy
    double _R;

    //! Exciton effective mass
    double _m;

    //! ISC time
    double _t_isc;

    //! Singlet diffusion coefficient
    double _sD;

    //! Triplet diffusion coefficient
    double _tD;

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
ExcitonSimple::~ExcitonSimple(void)
{
}


inline
PhysicalModelInterface*
ExcitonSimple::create_new(void) const
{
  return new ExcitonSimple(get_options());
}


inline
ExcitonSimple*
ExcitonSimple::create(const ModelOptions& options)
{
  return new ExcitonSimple(options);
}


#endif /* _ExcitonSimple_H_*/
