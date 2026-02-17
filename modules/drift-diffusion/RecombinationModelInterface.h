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
 * \file RecombinationModelInterface.h
 * \brief tiberCAD driftdiffusion module header.
 *
 * \note This file is part of module driftdiffusion.
 */


#ifndef TC_RECOMBINATIONMODELINTERFACE_H
#define TC_RECOMBINATIONMODELINTERFACE_H



#include "DriftDiffusionModelInterface.h"

#include <vector>

class Boundary;


//! Base class for recombination models
/*!
 * This is the base class for recombination models. A new recombination model
 * can be implemented by deriving from this class.
 */
class TBDLEXPORT RecombinationModelInterface : public DriftDiffusionModelInterface
{

  public:

    //! Destructor
    virtual ~RecombinationModelInterface(void);

    //! Get the electron and hole net recombination rates
    /*!
     * The net recombination rate is defined as \f$R_{net}=R-G\f$
     */
    virtual void get_net_recombination_rates(double& recomb_e,
        double& recomb_h) = 0;

    //! Get the electron and hole net recombination rate derivatives
    /*!
     * The net recombination rate is defined as \f$R_{net}=R-G\f$
     * The derivative is intended with respect to the electron (first
     * component) and hole (second component) density.
     */
    virtual void get_net_recombination_rate_derivatives(
        std::vector<double>& recomb_e, std::vector<double>& recomb_h) = 0;


    //! Creates a new named recombination model
    /*!
     * The model is created according to the given model name.
     * If it is not known, the NULL pointer is returned.
     *
     * \param name the model name
     * \param options the options as given in the input file
     * \return a pointer to the newly created object
     */
    static RecombinationModelInterface* create(const std::string& name,
        const PhysicalObject* owner, const ModelOptions& options = ModelOptions());


    //! Get the associated tunneling contact pointer, or NULL
    const Boundary* get_tunneling_contact(void);


  protected:

    //! \copydoc DriftDiffusionProperties::DriftDiffusionProperties()
    RecombinationModelInterface(const ModelOptions& options);

    //! Set the associated contact for tunneling
    /*!
     * This method should be called if the recombination model
     * is modelling carrier tunneling to or from a contact.
     */
    void set_tunneling_contact(const Boundary* bd);


  private:

    //! The associated tunneling contact, or NULL
    const Boundary* _tunneling_boundary;

};


//
// inline methods
//

inline
RecombinationModelInterface::RecombinationModelInterface(const ModelOptions& options)
 : DriftDiffusionModelInterface(options),
   _tunneling_boundary(NULL)
{
}

inline
RecombinationModelInterface::~RecombinationModelInterface(void)
{
}

inline
const Boundary*
RecombinationModelInterface::get_tunneling_contact(void)
{
  return _tunneling_boundary;
}

inline
void
RecombinationModelInterface::set_tunneling_contact(const Boundary* bd)
{
  _tunneling_boundary = bd;
}



inline
RecombinationModelInterface*
RecombinationModelInterface::create(const std::string& name,
    const PhysicalObject* owner, const ModelOptions& options)
{
  return dynamic_cast<RecombinationModelInterface*>(
      PhysicalModel::create("recombination_" + name, owner, options));
}



#endif // TC_RECOMBINATIONMODELINTERFACE_H
