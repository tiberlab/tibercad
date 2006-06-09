// $Id$

#ifndef _EXCITONMODEL_H_
#define _EXCITONMODEL_H_


#include "SimulationOptions.h"
#include "ExcitonProperties.h"
#include "TiberCad.h"

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

    /*! \copydoc ExcitonProperties::calculate_all()
     * 
     * This implementation models the most simple exciton model
     */
    virtual void calculate_all(double fermi_x, const Point& coord);

    virtual void read_database(const Dummy&);

    //! Set the exciton recombination time
    void set_recombination_time(double tau)
      { _t = tau; };

    //! Set the exciton binding energy
    void set_binding_energy(double R)
      { _R = R; };

    //! Set the exciton effective mass
    void set_effective_mass(double m)
      { _m = m; };

    //! Set the exciton effective mass
    void set_mobility(double mu)
      { _mu = mu; };
      
  protected:

    //! \copydoc ExcitonProperties::prepare_element_data()
    virtual void prepare_element_data(void);

  private:

    //! Exciton recombination time
    double _t;

    //! Exciton binding energy
    double _R;

    //! Exciton effective mass
    double _m;

    //! Exciton mobility
    double _mu;

    //! The effective density of states
    double _DOS;


};


//
// inline members
//

inline
ExcitonModel::~ExcitonModel(void)
{
}




#endif /* _EXCITONMODEL_H_*/
