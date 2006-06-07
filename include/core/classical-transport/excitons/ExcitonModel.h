// $Id$

#ifndef _EXCITONMODEL_H_
#define _EXCITONMODEL_H_

#include "vector_value.h"

#include "SimulationOptions.h"
#include "PhysicalProperties.h"
#include "DriftDiffusionDefs.h"
#include "TiberCad.h"

// GNU scientific library
//#include <gsl/gsl_sf_fermi_dirac.h>

#include <vector>

// forward declarations
class Point;
class Elem;

class ExcitonModel : public ExcitonProperties
{
    
  public:
  
    //! The empty constructor.
    ExcitonModel(DriftDiffusion* dd);
       
    //! A default (empty) destructor.
    virtual ~ExcitonModel(void) {};

    /*! \copydoc ExcitonProperties::calculate_all()
     * 
     * This implementation models the most simple exciton model
     */
    virtual void calculate_all(double fermi_x, const Point& coord);

    //! Set the exciton recombination time
    void set_recombination_time(double tau)
      { _t = tau; };

    //! Set the exciton binding energy
    void set_binding_energy(double R)
      { _R = R; };
      
  protected:

    //! \copydoc ExcitonProperties::prepare_element_data()
    virtual void prepare_element_data(void);

  private:

    //! The Drift-Diffusion object to get band data and densities from
    DriftDiffusion* _drift_diffusion;

    //! Exciton recombination time
    double _t;

    //! Exciton binding energy;
    double _R;

    /**!
     * The factor 3 * pow(2 * PI / h^2)^1.5
     * for calculating the effective density of states
     */
    double _DOS_factor;


};


//
// inline members
//

inline
ExcitonModel::~ExcitonModel(void)
{
}




#endif /* _EXCITONMODEL_H_*/
