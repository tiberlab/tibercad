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
 * \file QuantumDOS.h
 * \brief tiberCAD common module header.
 *
 * \note This file is part of module common.
 */


#ifndef _QUANTUMDOS_H_
#define _QUANTUMDOS_H_


#include "tibercad/physics/misc/DensityOfStates.h"
#include "tibercad/base/IDSet.h"

class SimulationInterface;
class Embracing;

/*!
 * \brief Base class for density of states
 */
class TBDLEXPORT QuantumDOS : public DensityOfStates
{

  public:

    //! Destructor
    virtual ~QuantumDOS(void) {};


    //! Creator function
    static QuantumDOS* create(const ModelOptions& options);


    //! Add the name of a quantum density calculation for this particle
    /*!
     * Instead of calculating the statistical classical densities one can
     * use a quantum density calculation. There can be more than one
     * quantum density calculation for the same particle as one could consider
     * eg. different valleys using different models.
     *
     * \param name the name of the quantum density calculation
     */
    void add_quantum_density(const std::string& name);


    //! Get a pointer to the quantum density simulation
    SimulationInterface* get_quantum_simulation(void);


    //! Set up an embracing region
    void set_embracing(Embracing* embracing);



  protected:

    //! Constructor
    QuantumDOS(const ModelOptions& options);

    virtual void do_init(void);

    virtual void do_reinit(const Elem* elem);

    virtual void prepare_submodels(void);

    virtual void do_print_info(void);


    //! Get occupied states and derivative w.r.t. Ef
    virtual void
    calculate_density_and_derivative(std::vector<double>& result, double Ef, double Epot,
        double kT, double kTlattice, const Elem* elem, const Point& p) const;

    //overloading for Trap.C
    virtual void
    calculate_density_and_derivative(std::vector<double>& result, double Ef, double Epot, double kT, double kTlattice) const;

  private:


    //! The quantum density calculation
    std::vector<SimulationInterface*> _quantum_density;


    //! The region IDs of the barrier
    IDSet _barrier_ids;


    //! The ID of the Density variable
    std::vector<ID> _density_ids;


    //! The ID for the BandEdge3D variable
    std::vector<ID> _3D_edge;


    //! If \c true a continuum will be added on top of the quantum density
    bool _add_continuum;


    //! The embracing of classical and quantum calculation
    Embracing* _embracing;


    //! The classical DOS, where QDOS is not defined
    DensityOfStates* _classical;


};

//
// inline methods
//

inline
QuantumDOS*
QuantumDOS::create(const ModelOptions& options)
{
  return new QuantumDOS(options);
}


inline
SimulationInterface*
QuantumDOS::get_quantum_simulation(void)
{
  if (_quantum_density.size() > 0)
    return _quantum_density[0];

  return NULL;
}



#endif // _DELTADOS_H_
