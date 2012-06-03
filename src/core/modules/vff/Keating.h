#ifndef _KEATING_H_
#define _KEATING_H_

#include "PhysicalModelInterface.h"
#include "Material.h"

//! The base class for Keating model parameters
class Keating : public PhysicalModelInterface
{
  public:

  virtual ~Keating(void) {};

  //! Get alpha parameter
  const double get_alpha_0() const;

  //! Get the "c-direction" alpha parameter in case of wz
  const double get_alpha_1() const;

  //! Get the beta parameter
  const double get_beta_0() const;

  //! Get the "c-direction" beta parameter in case of wz
  const double get_beta_1() const;


  protected:

  //! Constructor
  Keating(const ModelOptions& options);

  double _alpha_0;
  double _beta_0;
  double _alpha_1;
  double _beta_1;

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
Keating::Keating(const ModelOptions& options) :
_alpha_0(0.0),
_alpha_1(0.0),
_beta_0(0.0),
_beta_1(0.0),
  PhysicalModelInterface(options)
{
}

#endif
