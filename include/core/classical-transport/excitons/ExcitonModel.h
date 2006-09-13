// $Id$

#ifndef _EXCITONMODEL_H_
#define _EXCITONMODEL_H_

#include "SimulationOptions.h"
#include "ExcitonProperties.h"
#include "TiberCad.h"
#include "TypeDefs.h"


#include <string>

// forward declarations
class DriftDiffusionProperties;
class Dummy;

class ExcitonModel : public ExcitonProperties
{
    
  public:
  
    //! The empty constructor.
    ExcitonModel(void);
       
    //! A default (empty) destructor.
    virtual ~ExcitonModel(void);

    //! \copydoc ExcitonProperties::calculate_densities()
    virtual void calculate_densities(double fermi_x);

    //! \copydoc ExcitonProperties::calculate_recombination_rate()
    virtual void calculate_recombination_rate(void);

    /*! \copydoc ExcitonProperties::calculate_all()
     * 
     * This implementation models the most simple exciton model
     */
    virtual void calculate_all(double fermi_x, const Point& coord);

    virtual double get_nonradiative_recombination_rate(void);

    virtual void read_database(const Dummy&);

    //! Set the model of exciton generation to be used
    void set_exciton_generation_model(const std::string& model_name);

    //! Set the exciton recombination times
    void set_recombination_times(double tau_r, double tau_nr = 1e100)
      { _t_r = tau_r; _t_nr = tau_nr; };

    //! Set the exciton binding energy
    void set_binding_energy(double R)
      { _R = R; };

    //! Set the exciton effective mass
    void set_effective_mass(double m)
      { _m = m; };

    //! Set the exciton mobility
    void set_mobility(double mu)
      { _mu = mu; };
      
  protected:

    //! \copydoc ExcitonProperties::prepare_element_data()
    virtual void prepare_element_data(void);

  private:

    //! Exciton recombination time
    double _t_r;
    double _t_nr;

    //! Exciton binding energy
    double _R;

    //! Exciton effective mass
    double _m;

    //! Exciton mobility
    double _mu;

    //! The effective density of states
    double _DOS;

    //! the ID of the exciton generation model
    ID _gen_mod_id;

    //! The name of the exciton generation model to be used
    std::string _exciton_generation_model;


};


//
// inline members
//

inline
ExcitonModel::~ExcitonModel(void)
{
}


inline
void
ExcitonModel::set_exciton_generation_model(const std::string& model_name)
{
  _exciton_generation_model = model_name;
}



#endif /* _EXCITONMODEL_H_*/
