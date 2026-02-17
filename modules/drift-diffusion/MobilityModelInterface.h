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
 * \file MobilityModelInterface.h
 * \brief tiberCAD driftdiffusion module header.
 *
 * \note This file is part of module driftdiffusion.
 */


#ifndef TC_MOBILITYMODELINTERFACE_H
#define TC_MOBILITYMODELINTERFACE_H



#include "DriftDiffusionModelInterface.h"

#include "vector_value.h"

#include <vector>

//! The base class for mobility models
class TBDLEXPORT MobilityModelInterface : public DriftDiffusionModelInterface
{

  public:

    //! Destructor
    virtual ~MobilityModelInterface(void);

    //! Get the mobility
    virtual double get_mobility(void) = 0;

    //! Get the derivative of the mobility w.r.t. the electric potential
    virtual double get_derivative_potential(void);

    //! Get the derivative of the mobility w.r.t. the gradient of the electric potential
    virtual void get_derivative_grad_potential(libMesh::RealGradient& dm);

    //! Get the derivatives with respect to the fermi-level gradient
    virtual void get_derivative_grad_fermi(libMesh::RealGradient& dm);

    //! Creates a new named mobility model
    /*!
     * The model is created according to the given model name.
     * If it is not known, the NULL pointer is returned.
     *
     * \param name the model name
     * \param options the options as given in the input file
     * \return a pointer to the newly created object
     */
    static MobilityModelInterface* create(const std::string& name,
        const Material* mat, const ModelOptions& options = ModelOptions());


    //! Set the carrier type
    /*!
     * \param type the carrier type, can be \c e or \c h
     */
    void set_carrier_type(char type);


    //! Get the carrier type
    char get_carrier_type(void) const;


  protected:

    //! Create a mobility model
    MobilityModelInterface(const ModelOptions& options);



  private:

    //! The type of carriers this model is for
    /*!
     * Can be \c e or \c h
     */
    char _carrier;

};



//
// inline methods
//


inline
MobilityModelInterface::MobilityModelInterface(const ModelOptions& options)
  : DriftDiffusionModelInterface(options)
{
  std::string p = get_options().get_option("particle", "electron");
  if (p == "electron")
    _carrier = 'e';
  else
    _carrier = 'h';
}


inline
MobilityModelInterface::~MobilityModelInterface(void)
{
}



inline
void
MobilityModelInterface::set_carrier_type(char type)
{
  _carrier = type;
}


inline
char
MobilityModelInterface::get_carrier_type(void) const
{
  return _carrier;
}

inline
double
MobilityModelInterface::get_derivative_potential(void)
{
  return 0.0;
}

inline
void
MobilityModelInterface::get_derivative_grad_potential(libMesh::RealGradient& dm)
{
  dm.zero();
}

inline
void
MobilityModelInterface::get_derivative_grad_fermi(libMesh::RealGradient& dm)
{
  dm.zero();
}





#endif // TC_MOBILITYMODELINTERFACE_H
