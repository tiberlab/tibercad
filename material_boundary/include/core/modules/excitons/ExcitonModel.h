// $Id$

#ifndef _EXCITONMODEL_H_
#define _EXCITONMODEL_H_

#include "SimulationOptions.h"
#include "ExcitonProperties.h"
#include "TiberCad.h"
#include "TypeDefs.h"


#include <string>
#include <map>

class SimulationInterface;

//! A simple Exciton model
class ExcitonModel : public ExcitonProperties
{
    
  public:
  
    //! The empty constructor.
    ExcitonModel(void);
       
    //! A default (empty) destructor.
    virtual ~ExcitonModel(void);

    //! This method creates an ExcitonModel object
    static ExcitonModel* create(void);


  protected:

    /*! \copydoc ExcitonProperties::prepare_element_data() */
    virtual void prepare_element_data(void);

    /*! \copydoc PhysicalModel::create_new() */
    virtual PhysicalModelInterface* create_new(void) const;

    /*! \copydoc PhysicalModel::do_init_alloy() */
    virtual void do_init_alloy(const PhysicalModelInterface* comp_A,
        const PhysicalModelInterface* comp_B, double xa);

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

  private:

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
  return new ExcitonModel();
}


inline
ExcitonModel*
ExcitonModel::create(void)
{
  return new ExcitonModel();
}




#endif /* _EXCITONMODEL_H_*/
