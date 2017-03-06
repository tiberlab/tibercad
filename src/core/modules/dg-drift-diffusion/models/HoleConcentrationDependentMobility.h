/*
 * HoleConcentrationDependentMobility.h
 *
 *  Created on: Oct 11, 2016
 *      Author: mpatria
 */

#ifndef HOLECONCENTRATIONDEPENDENTMOBILITY_H_
#define HOLECONCENTRATIONDEPENDENTMOBILITY_H_

#include "MobilityModelInterface.h"


//! Sheet hole concentration dependent mobility model

//! It's used to simulate the behavior of the hole mobility dependent by sheet hole concentration in
//! diamond device with h-terminated diamond layer under the contact

/* The formula has been obtained by fitting experimental data given us  by Alta Frequenza group
 *
 * mu = d + {a - d} over {1 + (N_hc/c)^b}
 *
 * The constant a,b,c,d are obtained by using symmetric sigmoid function to fit the point data, and their values
 * dependent by the curve shape and the material
 *
 */

class TBDLLOCAL HoleConcentrationDependentMobility : public MobilityModelInterface
{

  public:

    //! Destructor
    virtual ~HoleConcentrationDependentMobility(void);

    //! Create a HoleConcentrationDependentMobility object
    static HoleConcentrationDependentMobility* create(const ModelOptions& options);

    //! \copydoc MobilityModelInterface::get_mobility()
    virtual double get_mobility(void);

    //! \copydoc MobilityModelInterface::get_mobility_derivatives()
    virtual void get_mobility_derivatives(std::vector<double>& dm);


  protected:

    //! constructor
    HoleConcentrationDependentMobility(const ModelOptions& options);

    //! \copydoc MobilityModelInterface::read_database()
    virtual void read_database(void);

    //! \copydoc MobilityModelInterface::do_init()
    virtual void do_init(void);

    //! \copydoc MobilityModelInterface::create_new()
    virtual PhysicalModelInterface* create_new(void) const;


  private:

    double a_;

    double b_;

    double c_;

    double d_;

};

//
// inline methods
//

inline
HoleConcentrationDependentMobility::HoleConcentrationDependentMobility(const ModelOptions& options)
  : MobilityModelInterface(options),
    a_(1),
    b_(1),
    c_(1),
    d_(1)

{
}


inline
HoleConcentrationDependentMobility*
HoleConcentrationDependentMobility::create(const ModelOptions& options)
{
  return new HoleConcentrationDependentMobility(options);
}


inline
PhysicalModelInterface*
HoleConcentrationDependentMobility::create_new(void) const
{
  return new HoleConcentrationDependentMobility(get_options());
}


inline
HoleConcentrationDependentMobility::~HoleConcentrationDependentMobility(void)
{
}


#endif /* HOLECONCENTRATIONDEPENDENTMOBILITY_H_ */
