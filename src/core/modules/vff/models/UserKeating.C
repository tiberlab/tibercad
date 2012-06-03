#include "UserKeating.h"

TIBER_MODULE(UserKeating, keating, user)

UserKeating::UserKeating(const ModelOptions& options):
Keating(options)
{

}

void
UserKeating::do_init(void)
{
  assign_alpha();
  assign_beta();
}

void
UserKeating::assign_alpha(void)
{
  if (get_material()->get_structure() == "zb")
    {
      if (get_options().find_option("alpha"))
        {
          _alpha_0 = get_option("alpha", 0.0);
          std::cout << "set alpha " << _alpha_0;
          _alpha_1 = get_option("alpha", 0.0);
        }
      else
        {
          std::string msg("Badly defined Keating parameters.");
          msg += "One parameter alpha is needed for zincoblenda";
          throw InitFailedException(msg);
        }
    }
  else if (get_material()->get_structure() == "wz")
    {
      if (get_options().find_option("alpha"))
        {
          _alpha_0 = get_option("alpha", 0.0);
          _alpha_1 = get_option("alpha", 0.0);
        }
      else if (get_options().find_option("alpha_0") &&
          get_options().find_option("alpha_1"))
        {
          _alpha_0 = get_option("alpha_0", 0.0);
          _alpha_1 = get_option("alpha_1", 0.0);
        }
      else
        {
          std::string msg("Badly defined Keating parameters.");
          msg += "One parameter alpha or two parameter alpha_0, alpha_1 are needed for wurtzite";
          throw InitFailedException(msg);
        }

    }

}

void
UserKeating::assign_beta(void)
{
  if (get_material()->get_structure() == "zb")
    {
      if (get_options().find_option("beta"))
        {
          _beta_0 = get_option("beta", 0.0);
          _beta_1 = get_option("beta", 0.0);
        }
      else
        {
          std::string msg("Badly defined Keating parameters.");
          msg += "One parameter beta is needed for zincoblenda";
          throw InitFailedException(msg);
        }
    }
  else if (get_material()->get_structure() == "wz")
    {
      if (get_options().find_option("beta"))
        {
          _beta_0 = get_option("beta", 0.0);
          _beta_1 = get_option("beta", 0.0);
        }
      else if (get_options().find_option("beta_0") &&
          get_options().find_option("beta_1"))
        {
          _beta_0 = get_option("beta_0", 0.0);
          _beta_1 = get_option("beta_1", 0.0);
        }
      else
        {
          std::string msg("Badly defined Keating parameters.");
          msg += "One parameter beta or two parameter beta_0, beta_1 are needed for wurtzite";
          throw InitFailedException(msg);
        }

    }

}
