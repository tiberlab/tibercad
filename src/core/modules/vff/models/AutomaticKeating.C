#include "AutomaticKeating.h"
#include "Database.h"
#include "Messages.h"
#include "RuntimeException.h"

TIBER_MODULE(AutomaticKeating, keating, automatic)

AutomaticKeating::AutomaticKeating(const ModelOptions& options):
Keating(options),
_c11(0.0),
_c12(0.0),
_c44(0.0),
_c13(0.0),
_c33(0.0)
{

}

void
AutomaticKeating::do_init(void)
{
  Keating::do_init();

  std::string warning1("Calculating keating parameters for alloy with Vegard's law. This is not safe.");
  if (get_material()->get_structure() == "zb")
    {
      parse_zb_database();
      if ((_alpha_0 == 0.0) || (_alpha_1 == 0.0))
        {
          if (get_material()->is_alloy())
            Messages::warning(warning1);
          calculate_zb_alpha();
        }
      if ((_beta_0 == 0.0) || (_beta_1 == 0.0))
        {
          if (get_material()->is_alloy())
            Messages::warning(warning1);
          calculate_zb_alpha();
        }
    }

  if (get_material()->get_structure() == "wz")
    {
      parse_wz_database();
      if ((_alpha_0 == 0.0) || (_alpha_1 == 0.0))
        {
          if (get_material()->is_alloy())
            Messages::warning(warning1);
          calculate_wz_alpha();
        }
      if ((_beta_0 == 0.0) || (_beta_1 == 0.0))
        {
          if (get_material()->is_alloy())
            Messages::warning(warning1);
          calculate_wz_alpha();
        }
    }

}

void
AutomaticKeating::parse_zb_database(void)
{
  double alpha, beta;
  const Database& db = get_database();
  db.set_section("keating");

  alpha = db.get("alpha", 0.0, false);
  _alpha_0 = alpha;
  _alpha_1 = alpha;
  beta = db.get("beta", 0.0, false);
  _beta_0 = beta;
  _beta_1 = beta;

  db.set_section("elasticity");
  _c11 = db.get("C11", 0.0, true);
  _c12 = db.get("C12", 0.0, true);
  _c44 = db.get("C44", 0.0, true);


}


void
AutomaticKeating::parse_wz_database(void)
{
  double alpha, beta;
  const Database& db = get_database();
  db.set_section("keating");

  alpha = db.get("alpha", 0.0, false);
  _alpha_0 = alpha;
  _alpha_1 = alpha;
  beta = db.get("beta", 0.0, false);
  _beta_0 = beta;
  _beta_1 = beta;
  alpha = db.get("alpha_0", 0.0, false);
  _alpha_0 = alpha;
  alpha = db.get("alpha_1", 0.0, false);
  _alpha_1 = alpha;
  beta = db.get("beta_0", 0.0, false);
  _beta_0 = beta;
  beta = db.get("beta_1", 0.0, false);
  _beta_1 = beta;

  db.set_section("elasticity");
  _c11 = db.get("C11", 0.0, true);
  _c12 = db.get("C12", 0.0, true);
  _c44 = db.get("C44", 0.0, true);
  _c13 = db.get("C13", 0.0, true);
  _c33 = db.get("C33", 0.0, true);


}


void
AutomaticKeating::calculate_zb_alpha(void)
{
  double alpha = (_c11 + 3.0 * _c12) * (_a  / 4.0);
  _alpha_0 = alpha;
  _alpha_1 = alpha;
}


void
AutomaticKeating::calculate_zb_beta(void)
{

  double beta = (_c11 - _c12) *  (_a  / 4.0);
  _beta_0 = beta;
  _beta_1 = beta;
}


void
AutomaticKeating::calculate_wz_alpha(void)
{
  std::string msg("Cannot calculate wz keating alpha. If you reached this point, you did something wrong.");
  throw RuntimeException(msg);
}



void
AutomaticKeating::calculate_wz_beta(void)
{
  std::string msg("Cannot calculate wz keating beta. If you reached this point, you did something wrong.");
  throw RuntimeException(msg);
}

