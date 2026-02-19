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
 * \file DegradationH2O.h
 * \brief tiberCAD pvmodule module header.
 *
 * \note This file is part of module pvmodule.
 */

#ifndef TC_DEGRADATIONH2O_H
#define TC_DEGRADATIONH2O_H

#include "DegradationModel.h"
#include "tibercad/module/SolutionProvider.h"

/*!
 * \brief An example for degradation due to relative humidity
 *
 * This class implements an equiv. circuit parameter dependency 
 * on relative humidity, calculated from a water ingress model.
 * The formulas have been obtained by fitting to data in
 * Bhatt et al., Organic Electronics 39 (2016) 258e266.
 * 
 * The fit functions are power laws, in particular:
 * 
 * \f{eqnarray*}
 *  \frac{R_s}{R_{s,0}} & = & 1 + \left(\frac{RH}{RH_0})^\gamma \\
 *  \frac{I_ph}{I_{ph,0}} & = & \frac{1}{1 + \left(\frac{RH}{RH_0})^\gamma} \\
 * \f}
 * 
 * Here \f$ RH \f$ is the relative humidity, which has to be provided by
 * another module.
 */
class DegradationH2O : public DegradationModel
{

  public:

    virtual ~DegradationH2O(void) = default;

    static DegradationH2O* create(const ModelOptions& options);


  protected:

    //! Protected constructor
    DegradationH2O(const ModelOptions& options);

    virtual void do_init(void) override;

    virtual void do_degrade_params(const libMesh::Elem* elem,
                                   const libMesh::Point& p,
                                   DegradationModel::Parameters& params) const final;


  private:

    //! The reference humidity in the photocurrent degradation fit
    double _RH_ref_ph = 1e9;

    //! The exponent in the photocurrent degradation fit
    double _exponent_ph = 1;

    //! The reference humidity in the series resistance degradation fit
    double _RH_ref_rs = 1e9;

    //! The exponent in the series resistance degradation fit
    double _exponent_rs = 1;

    //! From where to get relative humidity
    SolutionProvider _humidity_model;

};


#endif // TC_DEGRADATIONH2O_H
