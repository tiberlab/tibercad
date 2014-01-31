// $Id$

#include "HamiltonianModel.h"

class EfaModel : public HamiltonianModel 
{

  public:

    virtual ~EfaModel(void){};

    static EfaModel* create(const ModelOptions& options);

  protected:

    EfaModel(const ModelOptions& options);

    void do_init(void);

};

inline
EfaModel*
EfaModel::create(const ModelOptions& options)
{
  return new EfaModel(options);
}


