/*  
 * This file is part of the tiberCAD module driftdiffusion.
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
 * \file SchottkyTunneling.h
 * \brief tiberCAD driftdiffusion module header.
 *
 * \note This file is part of module driftdiffusion.
 */


#ifndef TC_SCHOTTKYTUNNELING_H
#define TC_SCHOTTKYTUNNELING_H

#include "RecombinationModelInterface.h"
#include "tibercad/base/HashMap.h"
#include "tibercad/base/TypeDefs.h"

#include "elem.h"
#include "point.h"

class SimulationInterface;
//class Elem;

//! Implementation of local Schottky Tunneling model
/*!
 * This class implements Schottky tunneling casted into
 * a local recombination model
 */
class TC_DLLOCAL SchottkyTunneling : public RecombinationModelInterface
{

  public:

    //! Destructor
    virtual ~SchottkyTunneling(void) {};

    //! Create a ConstantMobility object
    static SchottkyTunneling* create(const ModelOptions& options);

    //! \copydoc RecombinationModelInterface::get_net_recombination_rates()
    void get_net_recombination_rates(double& recomb_e, double& recomb_h);

    /*!
     * \copydoc
     * RecombinationModelInterface::get_net_recombination_rate_derivatives()
     */
    void get_net_recombination_rate_derivatives(
        std::vector<double>& recomb_e, std::vector<double>& recomb_h);


    
  protected:

    //! Constructor
    SchottkyTunneling(const ModelOptions& options);

    //! \copydoc RecombinationModelInterface::read_database()
    virtual void read_database(void);

    //! \copydoc RecombinationModelInterface::do_init()
    virtual void do_init(void);

    // ! \copydoc RecombinationModelInterface::do_reinit()
    //virtual void do_reinit(void);


  private:


    //! The maximum tunneling length in nm
    double _max_tunnel_length;


    //! The name of the associated contact
    std::string _contact_name;


    //! The contact voltage
    double _contact_voltage;


    //! The band
    char _band;


    //! The Schottky barrier
    double _barrier;


    //! The effective mass
    double _mass;


    //! A map containing all elements inside the tunnel length
    HashMap<const libMesh::Elem*, libMesh::Point>::Type _elem_map;



};



//
// inline methods
// 



inline
SchottkyTunneling*
SchottkyTunneling::create(const ModelOptions& options)
{
  return new SchottkyTunneling(options);
}






#endif // TC_SCHOTTKYTUNNELING_H
