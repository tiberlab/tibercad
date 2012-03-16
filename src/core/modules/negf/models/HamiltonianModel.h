// $Id$

#include "PhysicalModelInterface.h"
#include "tensor_value.h"
#include "vector_value.h"

class HamiltonianModel : public PhysicalModelInterface
{

  public:

    virtual ~HamiltonianModel(void) {};

    static HamiltonianModel* create(const ModelOptions& options);

    const TensorValue<double>& get_inv_mass(void) const;

    double get_degeneracy(void) const;

    const std::string& get_band_type(void) const;

  protected:

    HamiltonianModel(const ModelOptions& options);

    virtual void read_database(void);

    virtual void do_init(void);

    void set_invmass_tensor(void);

  private:

    TensorValue<double> _inv_mass_crys;

    TensorValue<double> _inv_mass;

    double _degeneracy;

    std::string _band_type;

};

inline
HamiltonianModel*
HamiltonianModel::create(const ModelOptions& options)
{
  return new HamiltonianModel(options);
}

inline
const TensorValue<double>&
HamiltonianModel::get_inv_mass(void) const
{
  return _inv_mass;
}

inline
double
HamiltonianModel::get_degeneracy(void) const
{
  return _degeneracy;
}

inline
const std::string&
HamiltonianModel::get_band_type(void) const
{
  return _band_type;
}

