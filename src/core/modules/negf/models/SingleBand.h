// $Id$

#include "PhysicalModelInterface.h"
#include "HamiltonianModel.h"
#include "tensor_value.h"
#include "vector_value.h"

class SingleBand : public HamiltonianModel 
{

  public:

    virtual ~SingleBand(void) {};

    static SingleBand* create(const ModelOptions& options);

  protected:

    SingleBand(const ModelOptions& options);

    void read_database(void);

    void do_init(void);

    void set_invmass_tensor(void);

};

inline
SingleBand*
SingleBand::create(const ModelOptions& options)
{
  return new SingleBand(options);
}


