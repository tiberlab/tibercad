#ifndef _KEATING_H_
#define _KEATING_H_

#include "PhysicalModel.h"
#include "Material.h"
#include "Messages.h"

//! The base class for Keating model parameters
class Keating : public PhysicalModel
{
  public:

    virtual ~Keating(void) {};

    //! Constructor
    Keating(const ModelOptions& options);

    //! Get alpha parameter
    double get_alpha_0() const;

    //! Get the "c-direction" alpha parameter in case of wz
    double get_alpha_1() const;

    //! Get the beta parameter
    double get_beta_0() const;

    //! Get the "c-direction" beta parameter in case of wz
    double get_beta_1() const;

    //! Get the distance parameter
    double get_d_0() const;

    //! Get the "c-direction" distance parameter in case of wz
    double get_d_1() const;

    //! Get the angle parameter
    double get_costeta_0() const;

    //! Get the "c-direction" angle parameter in case of wz
    double get_costeta_1() const;

  protected:

    //! Init operation common to all derived classes
    void do_init(void);


    double get_a(void) { return _a; }
    double get_c(void) { return _c; }
    double get_u(void) { return _u; }

    double& alpha_0(void) { return _alpha_0; }
    double& alpha_1(void) { return _alpha_1; }
    double& beta_0(void) { return _beta_0; }
    double& beta_1(void) { return _beta_1; }

    virtual void do_print_info(void);

  private:

    double _d_0;
    double _d_1;
    double _costeta_0;
    double _costeta_1;

    double _alpha_0;
    double _beta_0;
    double _alpha_1;
    double _beta_1;

    double _a;
    double _c;
    double _u;

    std::string _structure;


};

inline
double
Keating::get_alpha_0(void) const
{
  //std::string msg("Keating parameters alpha is 0");
  //if ((_alpha_0 == 0.0))
  //  Messages::warning(msg);
  return _alpha_0;
}

inline
double
Keating::get_alpha_1(void) const
{
  //std::string msg("Keating parameters alpha is 0");
  //if ((_alpha_1 == 0.0))
  //  Messages::warning(msg);
  return _alpha_1;
}

inline
double
Keating::get_beta_0(void) const
{
  std::string msg("Keating parameters beta is 0");
//  std::cout << "beta " << _beta_0 << "material " << get_material()->get_name();
  //if ((_beta_0 == 0.0))
  //  Messages::warning(msg);
  return _beta_0;
}

inline
double
Keating::get_beta_1(void) const
{
  std::string msg("Keating parameters beta is 0");
//  std::cout << "beta " << _beta_0 << "material " << get_material()->get_name();
  //if ((_beta_1 == 0.0))
  //  Messages::warning(msg);
  return _beta_1;
}

inline
double
Keating::get_d_0(void) const
{
  std::string msg("Vff parameters d is 0");
  //if ((_d_0 == 0.0))
  //  Messages::warning(msg);
  return _d_0;
}

inline
double
Keating::get_d_1(void) const
{
  std::string msg("Vff parameters d is 0");
  //if ((_d_1 == 0.0))
  //  Messages::warning(msg);
  return _d_1;
}

inline
double
Keating::get_costeta_0(void) const
{
  std::string msg("Vff parameters teta is 0");
  //if ((_costeta_0 == 0.0))
  //  Messages::warning(msg);
  return _costeta_0;
}

inline
double
Keating::get_costeta_1(void) const
{
  std::string msg("Vff parameters teta is 0");
  //if ((_costeta_1 == 0.0))
  //  Messages::warning(msg);
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
PhysicalModel(options)
{
}

#endif
