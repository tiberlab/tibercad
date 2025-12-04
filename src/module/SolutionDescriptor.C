/*
 * This file is part of tiberCAD.
 *
 * tiberCAD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * tiberCAD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with tiberCAD. If not, see <https://www.gnu.org/licenses/>.
 */

/*!
 * \file SolutionDescriptor.C
 * \brief tiberCAD API implementation.
 */


#include "tibercad/module/SolutionDescriptor.h"
#include "tibercad/utils/Utils.h"

#define CASE(key, str) case SolutionDescriptor::key: os << #str; break;
#define DEFAULT default: os << "unknown"; break;


std::ostream& operator<<(std::ostream& os, SolutionDescriptor::Type type)
{
  switch (type)
  {
    CASE(REAL, real)
    CASE(COMPLEX, complex)
    CASE(VECTOR, vector)
    CASE(TENSOR, tensor)
    CASE(NTUPLE, n-tuple)
    DEFAULT
  }
  return os;
}


std::ostream& operator<<(std::ostream& os, SolutionDescriptor::Location location)
{
  switch (location)
  {
    CASE(NODES, nodal)
    CASE(CELL, cell)
    CASE(ATOM, atom)
    CASE(GLOBAL, none)
    DEFAULT
  }
  return os;
}


void
SolutionDescriptor::get_components(const std::string& str, std::set<int> comp) const
{
  std::vector<std::string> tokens;
  Utils::tokenize(str, tokens, ":");

  for (int i = 0; i < tokens.size(); ++i)
  {
    switch (_type)
    {
      case COMPLEX:
        if (tokens[i] == "real")
          comp.insert(0);
        else if (tokens[i] == "imag")
          comp.insert(1);
        break;

      case VECTOR:
        if      (tokens[i] == "x")
          comp.insert(0);
        else if (tokens[i] == "y")
          comp.insert(1);
        else if (tokens[i] == "z")
          comp.insert(2);
        break;

      case TENSOR:
        if      (tokens[i] == "xx")
          comp.insert(0);
        else if (tokens[i] == "yy")
          comp.insert(1);
        else if (tokens[i] == "zz")
          comp.insert(2);
        else if (tokens[i] == "xy")
          comp.insert(3);
        else if (tokens[i] == "yz")
          comp.insert(4);
        else if (tokens[i] == "xz")
          comp.insert(5);
        break;

      case NTUPLE:
      {
        int n = atoi(tokens[i].c_str());
        if (n < n_components())
          comp.insert(n);
        break;
      }

      default:
        break;
    }
  }


  if (comp.empty())
  {
    for (int i = 0; i < n_components(); ++i)
      comp.insert(i);
  }
}
