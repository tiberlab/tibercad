#ifndef  _ATOMISTICGENERATOR2D_H_
#define _ATOMISTICGENERATOR2D_H_

//---------------------------------------------------------------------------------------------


#include "AtomisticGenerator.h"


class AtomisticGenerator2D : public AtomisticGenerator
{


public:

 AtomisticGenerator2D(AtomisticStructure* const as);

 ~AtomisticGenerator2D(void);

  static AtomisticGenerator2D* create(AtomisticStructure* const as);

protected:

  virtual void build();

  //virtual void passivate();

private:

};














#endif
