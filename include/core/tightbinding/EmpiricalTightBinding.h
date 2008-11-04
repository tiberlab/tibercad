#ifndef _EMPIRICALTIGHTBINDING_H_
#define _EMPIRICALTIGHTBINIDNG_H_

#include "tiber_config.h"

#include "TightBinding.h"


class EmpiricalTightBinding : public TightBinding
{

public:

  static EmpiricalTightBinding* create();

};

// Inline members definition
//----------------------------

inline
EmpiricalTightBinding* EmpiricalTightBinding::create()
{
  return new EmpiricalTightBinding;
}










#endif
