#ifndef  _ATOMISTICGENERATOR1D_H_
#define _ATOMISTICGENERATOR1D_H_

//---------------------------------------------------------------------------------------------

#include "AtomisticGenerator.h"


class AtomisticGenerator1D : public AtomisticGenerator 
{


public:

 AtomisticGenerator1D(AtomisticStructure* const as);

 ~AtomisticGenerator1D(void);

  static AtomisticGenerator1D* create(AtomisticStructure* const as);

protected:

  virtual void build();

  virtual void passivate();

private:

};









#endif // _ATOMISTICGENERATOR1D_H_
