#ifndef _STRAINEDSEMICONDUCTOR_H_
#define _STRAINEDSEMICONDUCTOR_H_

#include "SemiconductorModel.h"

class Elem;
class Macrostrain;

//! A drift-diffusion model for a strained semiconductor
/*!
 * This implementation of a strained semiconductor is based on
 * the \c SemiconductorModel class. Based on a \c Macrostrain object
 * it will recalculate for every element the most important equilibrium
 * properties and consider strain induced polarization.
 */
class StrainedSemiconductorModel : public SemiconductorModel
{

  public:
    StrainedSemiconductorModel(Macrostrain* strain);
    StrainedSemiconductorModel(const StrainedSemiconductorModel& model);
    virtual ~StrainedSemiconductorModel(void) {};

    //! Ignore strain related effects
    void ignore_strain(void);

    //! Include strain related effects
    void include_strain(void);

    //! Calculate the equilibrium properties
    /*!
     * This method calculates the equilibrium material properties
     * (eq. quasi Fermi level, eq. densities)
     */
    virtual void calculate_equilibrium_properties(
        int coupling = DriftDiffusionDefs::BOTH,
        double temperature = SimulationOptions::T);

  protected:

    virtual void prepare_element_data(void);

  private:

    Macrostrain* _strain;
    bool _ignore_strain;

    //! Factor to calculate strain induced correction to equilibrium
    //! Fermi level
    /*!
     * TODO describe mathematics
     */
    double gamma;

    //! The unstrained equilibrium band gap
    double Eg0;
    
    //! The unstrained equilibrium fermi level
    double Ef0;

    //! The unstrained equilibrium conduction band edge
    double Ec0;


};

inline
void
StrainedSemiconductorModel::ignore_strain(void)
{
  _ignore_strain = true;
}

inline
void
StrainedSemiconductorModel::include_strain(void)
{
  _ignore_strain = false;
}


#endif
