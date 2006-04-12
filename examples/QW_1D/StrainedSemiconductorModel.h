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

    void ignore_strain(void);

  protected:

    virtual void prepare_element_data(void);

  private:

    Macrostrain* _strain;
    bool _ignore_strain;


};

inline
void
StrainedSemiconductorModel::ignore_strain(void)
{
  _ignore_strain = true;
}

#endif
