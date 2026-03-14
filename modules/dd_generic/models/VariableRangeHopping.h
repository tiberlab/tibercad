/*  
 * This file is part of the tiberCAD module dd_generic.
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
 * \file VariableRangeHopping.h
 * \brief tiberCAD dd_generic module header.
 *
 * \note This file is part of module dd_generic.
 */

/*
 * VariableRangeHopping.h
 *
 *  Created on: 25 Jan 2022
 *      Author: miesu
 */

#ifndef TC_VARIABLERANGEHOPPING_H
#define TC_VARIABLERANGEHOPPING_H

#include "MobilityModelInterface.h"

 /*!
  * \brief Implementation of variable range hopping mobility
  *
  * Formula of the mobility
  *
  * $$
  * \mu=\mu_{0} \exp \left
  * \{-\left(\frac{k_{B} T_{h}}{k_{B} T+q E r}\right)^{\gamma}\right\}
  * $$
  *
  * Chen, T., van Gelder, J., van de Ven, B. et al. Classification
  * with a disordered dopant-atom network in silicon. Nature 577, 341–345 (2020).
  * https://doi.org/10.1038/s41586-019-1901-0
  */
class VariableRangeHopping : public MobilityModelInterface
{
  public:


    //! Destructor
    virtual ~VariableRangeHopping(void);


    //! \copydoc MobilityModelInterface::get_mobility()
    virtual double get_mobility(void) override;



    //! \copydoc MobilityModelInterface::get_derivative_grad_potential()
    /*!
     *   Derivative of the mobility
     *   $$
     *   \frac{d \mu}{d \nabla \phi}=-\frac{d \mu}{d \bar{E}}=-A
     *   \cdot B^{\gamma} \cdot \exp \left[-B^{\gamma}\right]
     *   \cdot \frac{\bar{E}}{|E|}
     *   $$
     *   \bigskip
     *   $$
     *   A=\frac{\mu_{0} \cdot q \cdot r \cdot \gamma}{k_{B}
     *   \cdot T_{h} + q \cdot |E| \cdot r}
     *   $$
     *   $$
     *   B=\frac{K_{B} T_{h}}{K_{B} T+q|E| \cdot r}
     *   $$
     */
    virtual void get_derivative_grad_potential(libMesh::RealGradient& dm);



  protected:



    //! Constructor
    VariableRangeHopping(const ModelOptions& options);



    //! \copydoc MobilityModelInterface::do_init()
    virtual void do_init(void) override;



  private:

    //! prefactor of mobility
    double _mu_0;

    //! average hopping distance r
    double _hopping_distance;

    //! characteristic temperature of hopping conduction Th
    double _temp_hopping_conduction;

    //! exponent
    double _gamma;
};









#endif /*_VARIABLERANGEHOPPING_H_ */
