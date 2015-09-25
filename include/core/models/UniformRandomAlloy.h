// $Id$

#ifndef _UNIFORMRANDOMALLOY_H_
#define _UNIFORMRANDOMALLOY_H_


#include "ExternalProfile.h"

#include <random>


//! Class to create a uniform random distribution
class UniformRandomAlloy : public ExternalProfile
{
  public:

    //! Constructor
    UniformRandomAlloy(const ModelOptions& options);

    //! Destructor
    virtual ~UniformRandomAlloy(void);

    virtual double get_data(const Elem* elem) const;

    virtual double get_data(const Elem* elem, const Point& p) const;

    virtual std::pair<double, double> get_min_max(void) const;

  private:


    double _min;
    double _max;

    double _mean_composition;

    //! The available site density in mesh units
    double _site_density;

    //! Random number generator engine
    mutable std::mt19937 _rnd_generator;

};


#endif // _UNIFORMRANDOMALLOY_H_
