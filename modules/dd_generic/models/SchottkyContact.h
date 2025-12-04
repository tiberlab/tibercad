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
 * \file SchottkyContact.h
 * \brief tiberCAD dd_generic module header.
 *
 * \note This file is part of module dd_generic.
 */


#ifndef _SCHOTTKYCONTACT_H_
#define _SCHOTTKYCONTACT_H_

#include "ElectricalContact.h"


/*!
 * \brief A Schottky contact
 */
class TBDLLOCAL SchottkyContact : public ElectricalContact
{
  public:

    //! Create a schottky contact
    static SchottkyContact* create(const ModelOptions& options);


  protected:

    //! The constructor
    SchottkyContact(const ModelOptions& options);

    /*! \copydoc ElectricalContact::do_init() */
    virtual void do_init(void);


    //! Calculate all coefficients
    virtual void do_compute(void);


  private:

    //! The metal Fermi level (negative of work function)
    double _metal_Ef = 0.0;

    //! The reference band
    std::string _band;

    //! Is this a fixed barrier or not?
    bool _fixed_barrier = true;


    //! Do we include thermionic emission?
    bool _thermionic_emission = true;

    /*!
     * \brief Do we use Scott and Malliaras field and mobility dependent recombination velocity?
     * For details, see Chem. Phys. Lett. 299 (1999) 115
     */
    bool _scott_malliaras = false;

    //! Image-force lowering
    bool _barrier_lowering = false;

    /*!
     * \brief Additional multiplicative factor
     * 
     * This can be used to adjust the Richardson constant, for example due to
     * additional thin oxide layers beneath the metal.
     */
    double _correction_factor = 1.0;

};



//
// inline
//


inline
SchottkyContact*
SchottkyContact::create(const ModelOptions& options)
{
  return new SchottkyContact(options);
}






#endif // _SCHOTTKYCONTACT_H_
