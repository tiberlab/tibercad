#include "UserKeating.h"
#include "PhysicalObject.h"

TIBER_MODULE(UserKeating, keating, user)

UserKeating::UserKeating(const ModelOptions& options):
Keating(options)
{

}

void
UserKeating::do_init(void)
{
  Keating::do_init();
  assign_alpha();
  assign_beta();
  assign_alpha_parents();
  assign_beta_parents();
  check_parameters();
}


void
UserKeating::assign_alpha_parents(void)
{
  ModelOptions::submodel_iterator it = get_options().submodels_begin("component");

  if (get_material()->get_structure() == "zb")
    {
      for (; it != get_options().submodels_end("component"); ++it)
        {
          if ((it->second).get_name() == get_material()->get_name()) // o simile
            {
              _alpha_0 = (it->second).get_option("alpha", 0.0);
              _alpha_1 = (it->second).get_option("alpha", 0.0);
            }
        }
    }

  if (get_material()->get_structure() == "wz")
    {
      for (; it != get_options().submodels_end("component"); ++it)
        {
          if ((it->second).get_name() == get_material()->get_name()) // o simile
            {
              _alpha_0 = (it->second).get_option("alpha", 0.0);
              _alpha_1 = (it->second).get_option("alpha", 0.0);
              _alpha_0 = (it->second).get_option("alpha_0", 0.0);
              _alpha_1 = (it->second).get_option("alpha_1", 0.0);
            }
        }
    }

}


void
UserKeating::assign_beta_parents(void)
{
  ModelOptions::submodel_iterator it = get_options().submodels_begin("component");

  if (get_material()->get_structure() == "zb")
    {
      for (; it != get_options().submodels_end("component"); ++it)
        {
          if ((it->second).get_name() == get_material()->get_name()) // o simile
            {
              _beta_0 = (it->second).get_option("beta", 0.0);
              _beta_1 = (it->second).get_option("beta", 0.0);
            }
        }
    }

  if (get_material()->get_structure() == "wz")
    {
      for (; it != get_options().submodels_end("component"); ++it)
        {
          if ((it->second).get_name() == get_material()->get_name()) // o simile
            {
              _beta_0 = (it->second).get_option("beta", 0.0);
              _beta_1 = (it->second).get_option("beta", 0.0);
              _beta_0 = (it->second).get_option("beta_0", 0.0);
              _beta_1 = (it->second).get_option("beta_1", 0.0);
            }
        }
    }

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


    }

}

void
UserKeating::check_parameters(void)
{
  std::string msg("Keating parameters are not correctly defined by user");
  if ((_alpha_0 == 0.0) || (_alpha_1 == 0.0) || (_beta_0 == 0.0) || (_beta_1 == 0.0))
    throw(msg);
}
