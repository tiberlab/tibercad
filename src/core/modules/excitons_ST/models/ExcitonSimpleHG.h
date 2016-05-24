
#ifndef _EXCITONSIMPLEHG_H_
#define _EXCITONSIMPLEHG_H_

#include "SimulationOptions.h"
#include "ExcitonProperties.h"
#include "TiberCad.h"
#include "TypeDefs.h"

#include <string>
#include <map>

class SimulationInterface;

//! A simple Exciton model with host-guest recombination
class ExcitonSimpleHG : public ExcitonProperties
{
    
  public:
  
    //! The empty constructor.
    ExcitonSimpleHG(const ModelOptions& options);
       
    //! A default (empty) destructor.
    virtual ~ExcitonSimpleHG(void);

    //! This method creates an ExcitonSimpleHG object
    static ExcitonSimpleHG* create(const ModelOptions& options);

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
    virtual double get_s_generation_rate(void);

    //! Get the singlet host-guest generation rate
    virtual double get_s_hg_generation_rate(void);

    //! Get the singlet host-guest recombination rate
    virtual double get_s_hg_recombination_rate(void);

    //! Get the triplet nonradiative recombination rate
    virtual double get_t_nonradiative_recombination_rate(void);

    //! Get the triplet radiative recombination rate
    virtual double get_t_radiative_recombination_rate(void);
 
    //! Get the triplet dissociation rate
    virtual double get_t_dissociation_rate(void);

    //! Get the triplet generation rate
    virtual double get_t_generation_rate(void);

    //! Get the triplet host-guest generation rate
    virtual double get_t_hg_generation_rate(void);

    //! Get the triplet host-guest recombination rate
    virtual double get_t_hg_recombination_rate(void);


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

    //! Material permittivity
    double _er;

    //! ISC time
    double _t_isc;

    //! Singlet diffusion coefficient
    double _sD;

    //! Triplet diffusion coefficient
    double _tD;

    //! Foster radius
    double _Rf;

    //! Dexter radius
    double _Rd;

    //! Host-guest average sites distance
    double _R_hg;

    //! Foster rate
    double _Kf;

    //! Dexter rate for singlets
    double _Kds;

    //! Dexter rate for triplets
    double _Kdt;

    //! The ID for the generation model to be used
    ID _gen_model;

    //! The ID for the exciton model to be used in host-guest simulations
    ID _exs_model;
    ID _ext_model;

    //! The ID for the band gap variable
    ID _Eg_id;

    //! The DriftDiffusion simulation to be used
    SimulationInterface* _dd_sim;

    //! The host-guest excitons simulation to be used
    SimulationInterface* _hg_sim;

    //! The band gap for each element
    std::map<const Elem*, double> _bandgap_data;


};


//
// inline members
//

inline
ExcitonSimpleHG::~ExcitonSimpleHG(void)
{
}


inline
PhysicalModelInterface*
ExcitonSimpleHG::create_new(void) const
{
  return new ExcitonSimpleHG(get_options());
}


inline
ExcitonSimpleHG*
ExcitonSimpleHG::create(const ModelOptions& options)
{
  return new ExcitonSimpleHG(options);
}


#endif
