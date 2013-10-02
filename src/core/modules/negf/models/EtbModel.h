// $Id$

#include "HamiltonianModel.h"

class EtbModel : public HamiltonianModel 
{

  public:

    virtual ~EtbModel(void){};

    static EtbModel* create(const ModelOptions& options);

  protected:

    EtbModel(const ModelOptions& options);

    void do_init(void);

};

inline
EtbModel*
EtbModel::create(const ModelOptions& options)
{
  return new EtbModel(options);
}


