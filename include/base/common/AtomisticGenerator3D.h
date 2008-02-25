#ifndef  _ATOMISTICGENERATOR3D_H_
#define _ATOMISTICGENERATOR3D_H_

//---------------------------------------------------------------------------------------------


#include "AtomisticGenerator.h"


class AtomisticGenerator3D : public AtomisticGenerator 
{


public:

 AtomisticGenerator3D(AtomisticStructure* const as);

 ~AtomisticGenerator3D(void);

  static AtomisticGenerator3D* create(AtomisticStructure* const as);

protected:

  virtual void build();

  virtual void passivate();

private:

};














#endif
