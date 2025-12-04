/*  
 * This file is part of the tiberCAD module pvmodule.
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
 * \file ElementaryPWL.C
 * \brief tiberCAD pvmodule module implementation.
 *
 * \note This file is part of module pvmodule.
 */

#include "ElementaryPWL.h"

#include "tibercad/module/TiberModule.h"

inline
ElementaryPWL::ElementaryPWL(const ModelOptions& options) :
  ElementaryCell(options)
{
}


inline
ElementaryPWL*
ElementaryPWL::create(const ModelOptions& options)
{
  ElementaryPWL* cd = new ElementaryPWL(options);

  return cd;
}



void
ElementaryPWL::do_init(void)
{
  // reading jv_ref.dat file, as reference for JV curve
  std::string jv_file;
  jv_file = get_option("jv_data", jv_file);
  if (jv_file == "")
    throw InitFailedException("JV of elementary cell needs to be "
                              "specified via 'jv_data' option for pwl elementary_cell model.");

  std::ifstream ifs;
  ifs.open(jv_file);
  if (ifs.fail() || !ifs.good())
    throw InitFailedException("Cannot read JV of elementary cell "
                              "from file " +
                              jv_file);

  size_t i = 0;
  const size_t buf_len = 256;
  char buf[buf_len];
  while (ifs.good())
  {
    if (i == _jv_v.size())
    {
      size_t n_new = _jv_v.size() + 10;
      _jv_v.reserve(n_new);
      _jv_j.reserve(n_new);
    }
    ifs.getline(buf, buf_len);
    if (buf[0] != '#')
    {
      std::stringstream in(buf);
      double l, s;
      if (in >> l >> s)
      {
        _jv_v.push_back(l);
        _jv_j.push_back(s);
        i++;
      }
    }
  }
  ifs.close();
}


void
ElementaryPWL::do_write_netlist(unsigned int top_node, unsigned int bottom_node,
                                unsigned int& /* next_free */,
                                double area,
                                const libMesh::Elem* /* elem */,
                                const libMesh::Point& /* p */,
                                std::ostream& os) const 
{
  os << "B" << top_node << " " << top_node << " " << bottom_node
     << " I=pwl(V(" << top_node << ")-V(" << bottom_node << ")";

  for (int n = 0; n < _jv_v.size(); n++)
    os << ", " << _jv_v[n] << ", " << _jv_j[n] * area << "m";

  os << ")\n\n";
}
