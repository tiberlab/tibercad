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
 * \file ExcitonISC.h
 * \brief tiberCAD dd_generic module header.
 *
 * \note This file is part of module dd_generic.
 */


#ifndef TC_EXCITONISC_H
#define TC_EXCITONISC_H

#include "RecombinationModelInterface.h"
#include "tibercad/base/TypeDefs.h"

class SimulationInterface;

//! Implementation of direct recombination
/*!
 * This class implements direct recombination processes that can be
 * modeled by \f[R_{direct}=C(np-n_i^2)\f]
 */
class TC_DLLOCAL ExcitonISC : public RecombinationModelInterface
{

  public:

    //! Destructor
    virtual ~ExcitonISC(void) {};

    //! Create a ConstantMobility object
    static ExcitonISC* create(const ModelOptions& options);


  protected:

    //! Constructor
    ExcitonISC(const ModelOptions& options);

    //! \copydoc RecombinationModelInterface::read_database()
    virtual void read_database(void) override;

    //! \copydoc RecombinationModelInterface::do_init()
    virtual void do_init(void) override;

    //! \copydoc RecombinationModelInterface::do_reinit()
    virtual void do_reinit(void) override;

    //! \copydoc RecombinationModelInterface::calculate_rate_and_derivatives()
    virtual void calculate_rate_and_derivatives(std::vector<double>& R, std::vector<std::vector<double>>& dPotentials) override;

  private:

    typedef std::map<std::pair<SimulationInterface*, SimulationInterface*>,
        std::pair<unsigned int, double> > QRecMap;

    //! Recombination rate parameter
    double  _C;

    //! The quantum optics simulation, if available
    SimulationInterface* _quantum_optics;

    //! The solution ID for the optical recombination
    ID _rec_id;

    //! A static map to put quantum recombination in
    /*!
     * This map is used so as to not calculate the same quantity
     * several times.
     */
    static QRecMap _qrec_vals;

};



//
// inline methods
// 

inline
ExcitonISC::ExcitonISC(const ModelOptions& options)
  : RecombinationModelInterface(options),
    _quantum_optics(NULL)
{
}


inline
ExcitonISC*
ExcitonISC::create(const ModelOptions& options)
{
  return new ExcitonISC(options);
}






#endif // TC_EXCITONISC_H
