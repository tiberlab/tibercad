/*  
 * This file is part of the tiberCAD module common.
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
 * \file BulkDOS.h
 * \brief tiberCAD common module header.
 *
 * \note This file is part of module common.
 */


#ifndef TC_BULKDOS_H
#define TC_BULKDOS_H


#include "tibercad/physics/misc/DensityOfStates.h"


/*!
 * \brief Base class for density of states
 */
class TBDLEXPORT BulkDOS : public DensityOfStates
{

  public:

    //! Destructor
    virtual ~BulkDOS(void) {};


    //! Creator function
    static BulkDOS* create(const ModelOptions& options);



  protected:

    //! Constructor
    BulkDOS(const ModelOptions& options);

    //! Read band edge, mass, degeneracy from database
    /*!
     * This will only read from the database if the \c particle option
     * is provided as \c electron or \c hole
     */
    virtual void read_database(void);

    virtual void do_init(void);

    virtual void do_print_info(void);


    //! Get occupied states and derivative w.r.t. Ef
    virtual void
    calculate_density_and_derivative(std::vector<double>& result, double Ef, double Epot,
        double kT, double kTlattice, const Elem* elem, const Point& p) const;


  private:

    //! A structure for band gap parameters
    struct GapParameters
    {
      double Ev;
      double Eg0;
      double varshni_a;
      double varshni_b;
      void zero(void) { Ev = Eg0 = varshni_a = varshni_b = 0.0; }
      double gap(double T) const { return Eg0 - varshni_a * T * T / (T + varshni_b); }
    };


    //! We allow for several subbands
    std::vector<double> _ref_energies;

    //! The effective masses
    std::vector<double> _dos_mass;

    //! The degeneracies
    std::vector<int> _degeneracy;

    //! The band gap parameters
    std::vector<GapParameters> _gap_params;


    //! The DOS factor
    /*!
     * In 3D the effective DOS is
     * \f[ N_{3D} = \f]
     */
    double _dos_factor;


};

//
// inline methods
//

inline
BulkDOS*
BulkDOS::create(const ModelOptions& options)
{
  return new BulkDOS(options);
}


#endif // TC_BULKDOS_H
