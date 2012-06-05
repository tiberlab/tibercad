#ifndef _KEATING_H_
#define _KEATING_H_

#include "PhysicalModelInterface.h"
#include "Material.h"
#include "RuntimeException.h"

//! The base class for Keating model parameters
class Keating : public PhysicalModelInterface
{
  public:

  virtual ~Keating(void) {};

  //! Constructor
  Keating(const ModelOptions& options);

  //! Get alpha parameter
  const double get_alpha_0() const;

  //! Get the "c-direction" alpha parameter in case of wz
  const double get_alpha_1() const;

  //! Get the beta parameter
  const double get_beta_0() const;

  //! Get the "c-direction" beta parameter in case of wz
  const double get_beta_1() const;

  //! Get the distance parameter
  const double get_d_0() const;

  //! Get the "c-direction" distance parameter in case of wz
  const double get_d_1() const;

  //! Get the angle parameter
  const double get_costeta_0() const;

  //! Get the "c-direction" angle parameter in case of wz
  const double get_costeta_1() const;

  protected:

  //! Init operation common to all derived classes
  void do_init(void);

  double _alpha_0;
  double _beta_0;
  double _alpha_1;
  double _beta_1;
  double _d_0;
  double _d_1;
  double _costeta_0;
  double _costeta_1;
  double _a;
  double _c;
  double _u;

  std::string _structure;


};

inline
const double
Keating::get_alpha_0(void) const
{
  std::string msg("Keating parameters alpha is 0");
  if ((_alpha_0 == 0.0))
    throw RuntimeException(msg);
  return _alpha_0;
}

inline
const double
Keating::get_alpha_1(void) const
{
  std::string msg("Keating parameters alpha is 0");
  if ((_alpha_1 == 0.0))
    throw RuntimeException(msg);
  return _alpha_1;
}

inline
const double
Keating::get_beta_0(void) const
{
  std::string msg("Keating parameters beta is 0");
//  std::cout << "beta " << _beta_0 << "material " << get_material()->get_name();
  if ((_beta_0 == 0.0))
    throw RuntimeException(msg);
  return _beta_0;
}

inline
const double
Keating::get_beta_1(void) const
{
  std::string msg("Keating parameters beta is 0");
//  std::cout << "beta " << _beta_0 << "material " << get_material()->get_name();
  if ((_beta_1 == 0.0))
    throw RuntimeException(msg);
  return _beta_1;
}

inline
const double
Keating::get_d_0(void) const
{
  std::string msg("Vff parameters d is 0");
  if ((_d_0 == 0.0))
    throw RuntimeException(msg);
  return _d_0;
}

inline
const double
Keating::get_d_1(void) const
{
  std::string msg("Vff parameters d is 0");
  if ((_d_1 == 0.0))
    throw RuntimeException(msg);
  return _d_1;
}

inline
const double
Keating::get_costeta_0(void) const
{
  std::string msg("Vff parameters teta is 0");
  if ((_costeta_0 == 0.0))
    throw RuntimeException(msg);
  return _costeta_0;
}

inline
const double
Keating::get_costeta_1(void) const
{
  std::string msg("Vff parameters teta is 0");
  if ((_costeta_1 == 0.0))
    throw RuntimeException(msg);
  return _costeta_1;
}

inline
Keating::Keating(const ModelOptions& options) :
_alpha_0(0.0),
_alpha_1(0.0),
_beta_0(0.0),
_beta_1(0.0),
_a(0.0),
_c(0.0),
_d_0(0.0),
_d_1(0.0),
_costeta_0(0.0),
_costeta_1(0.0),
_u(0.375),
PhysicalModelInterface(options)
{
}

#endif
