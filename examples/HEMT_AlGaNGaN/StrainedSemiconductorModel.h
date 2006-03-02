#ifndef _STRAINEDSEMICONDUCTOR_H_
#define _STRAINEDSEMICONDUCTOR_H_

#include "SimpleSemiconductorModel.h"

class Elem;
class Macrostrain;

class StrainedSemiconductorModel : public SimpleSemiconductorModel
{

  public:
    StrainedSemiconductorModel(Macrostrain* strain);
    StrainedSemiconductorModel(const StrainedSemiconductorModel& model);
    virtual ~StrainedSemiconductorModel(void) {};

    virtual void calculate_all(double potential, double Ef_e, double Ef_h,
        const Point& p, const Elem* elem,
        int coupling = DriftDiffusionDefs::BOTH);

  private:

    Macrostrain* _strain;

};

#endif
