#ifndef _KEATING_H_
#define _KEATING_H_

#include "PhysicalModelInterface.h"
#include "Material.h"

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
  const double get_teta_0() const;

  //! Get the "c-direction" angle parameter in case of wz
  const double get_teta_1() const;

  protected:

  //! Init operation common to all derived classes
  void do_init(void);

  double _alpha_0;
  double _beta_0;
  double _alpha_1;
  double _beta_1;
  double _d_0;
  double _d_1;
  double _teta_0;
  double _teta_1;
  double _a;
  double _c;
  double _u;

  std::string _structure;


};

inline
const double
Keating::get_alpha_0(void) const
{
  return _alpha_0;
}

inline
const double
Keating::get_alpha_1(void) const
{
  return _alpha_1;
}

inline
const double
Keating::get_beta_0(void) const
{
  return _beta_0;
}

inline
const double
Keating::get_beta_1(void) const
{
  return _beta_1;
}

inline
const double
Keating::get_d_0(void) const
{
  return _d_0;
}

inline
const double
Keating::get_d_1(void) const
{
  return _d_1;
}

inline
const double
Keating::get_teta_0(void) const
{
  return _teta_0;
}

inline
const double
Keating::get_teta_1(void) const
{
  return _teta_1;
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
_teta_0(0.0),
_teta_1(0.0),
_u(0.375),
PhysicalModelInterface(options)
{
}

#endif
