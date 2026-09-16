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
 * \file ElementaryPWL.h
 * \brief tiberCAD pvmodule module header.
 *
 * \note This file is part of module pvmodule.
 */

#ifndef TC_ELEMENTARYPWL_H
#define TC_ELEMENTARYPWL_H

#include "ElementaryCell.h"


class ElementaryPWL : public ElementaryCell
{

  public:

    virtual ~ElementaryPWL(void) override {};


  protected:

    ElementaryPWL(const ModelOptions& options);

    virtual void do_init(void) override;

    //! Write the netlist
    virtual void do_write_netlist(unsigned int top_node, unsigned int bottom_node,
                                  unsigned int& next_free,
                                  double area,
                                  const libMesh::Elem* elem,
                                  const libMesh::Point& p,
                                  std::ostream& os) const override;

  private:

    //! The voltage data points
    std::vector<double> _jv_v;

    //! The current density data points
    std::vector<double> _jv_j;

};


#endif // TC_ELEMENTARYPWL_H
