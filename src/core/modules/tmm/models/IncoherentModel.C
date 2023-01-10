
#include "IncoherentModel.h"
#include "TiberModule.h"


inline
IncoherentModel::IncoherentModel(const ModelOptions& options) :
  TmmBulkModel(options)
{
 _Incoheret_Index=0;
}

inline
IncoherentModel*
IncoherentModel::create(const ModelOptions& options)
{
  return new IncoherentModel(options);
}

inline
const double&
IncoherentModel::get_Incoherent_Index() const
{
 
  return _Incoheret_Index;
}

inline
void
IncoherentModel::set_Incoherent_Index(const double& Incoheret_Index)
{
  
  _Incoheret_Index =  Incoheret_Index;
  
}

