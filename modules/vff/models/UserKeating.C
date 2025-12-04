/*  
 * This file is part of the tiberCAD module vff.
 *
 * tiberCAD modules are licensed under the GNU General Public License v3.
 *
 * tiberCAD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published
 * by the Free Software Foundation, either version 3 of the License,
 * or (at your option) any later version.
 *
 * tiberCAD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with tiberCAD. If not, see <http://www.gnu.org/licenses/>.
 */

/*!
 * \file UserKeating.C
 * \brief tiberCAD vff module implementation.
 *
 * \note This file is part of module vff.
 */

#include "UserKeating.h"
#include "tibercad/physics/PhysicalObject.h"

#include "tibercad/module/TiberModule.h"


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
}


void
UserKeating::assign_alpha_parents(void)
{
  ModelOptions::submodel_iterator it = get_options().submodels_begin("component");

  if (get_material()->get_structure() == "zb")
    {
      for (; it != get_options().submodels_end("component"); ++it)
        {
          if ((it->second).get_name() == get_material()->get_name())
            {
              alpha_0() = (it->second).get_option("alpha", 0.0);
              alpha_1() = (it->second).get_option("alpha", 0.0);
            }
        }
    }

  if (get_material()->get_structure() == "wz")
    {
      for (; it != get_options().submodels_end("component"); ++it)
        {
          if ((it->second).get_name() == get_material()->get_name())
            {
              if ((it->second).find_option("alpha"))
                {
                  alpha_0() = (it->second).get_option("alpha", 0.0);
                  alpha_1() = (it->second).get_option("alpha", 0.0);
                }
              if ((it->second).find_option("alpha_0"))
                alpha_0() = (it->second).get_option("alpha_0", 0.0);
              if ((it->second).find_option("alpha_1"))
                alpha_1() = (it->second).get_option("alpha_1", 0.0);
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
              beta_0() = (it->second).get_option("beta", 0.0);
              beta_1() = (it->second).get_option("beta", 0.0);
            }
        }
    }

  if (get_material()->get_structure() == "wz")
    {
      for (; it != get_options().submodels_end("component"); ++it)
        {
          if ((it->second).get_name() == get_material()->get_name()) // o simile
            {
              if ((it->second).find_option("beta"))
                {
                  beta_0() = (it->second).get_option("beta", 0.0);
                  beta_1() = (it->second).get_option("beta", 0.0);
                }
              if ((it->second).find_option("beta_0"))
                beta_0() = (it->second).get_option("beta_0", 0.0);
              if ((it->second).find_option("beta_1"))
                beta_1() = (it->second).get_option("beta_1", 0.0);
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
          alpha_0() = get_option("alpha", 0.0);
          alpha_1() = get_option("alpha", 0.0);
        }

    }
  else if (get_material()->get_structure() == "wz")
    {
      if (get_options().find_option("alpha"))
        {
          alpha_0() = get_option("alpha", 0.0);
          alpha_1() = get_option("alpha", 0.0);
        }
      else if (get_options().find_option("alpha_0") &&
          get_options().find_option("alpha_1"))
        {
          alpha_0() = get_option("alpha_0", 0.0);
          alpha_1() = get_option("alpha_1", 0.0);
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
          beta_0() = get_option("beta", 0.0);
          beta_1() = get_option("beta", 0.0);
        }

    }
  else if (get_material()->get_structure() == "wz")
    {
      if (get_options().find_option("beta"))
        {
          beta_0() = get_option("beta", 0.0);
          beta_1() = get_option("beta", 0.0);
        }
      else if (get_options().find_option("beta_0") &&
          get_options().find_option("beta_1"))
        {
          beta_0() = get_option("beta_0", 0.0);
          beta_1() = get_option("beta_1", 0.0);
        }


    }

}
