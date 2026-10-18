/*  
 * This file is part of the tiberCAD module masstransport.
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
 * \file MassTransport.h
 * \brief tiberCAD masstransport module header.
 *
 * \note This file is part of module masstransport.
 */


#ifndef TC_MASSTRANSPORT_H
#define TC_MASSTRANSPORT_H

#include "tibercad/module/SimulationInterface.h"
#include "tibercad/solver/TiberLinearSystem.h"
#include <string>
using namespace std;


/*!
 * 
 * \brief This is a simple implementation of the Fick's laws for
 *        modeling mass transport e.g. into solar cells.
 *
 * The implementation uses the partial pressure as primary
 * variable, and assumes the solubility to be piecewise constant.
 */
class TC_DLLOCAL MassTransport : public SimulationInterface
{

  public:

    //! Destructor
    /*!
     * We do not declare it virtual here, as we will not allow
     * to derive from this class anyway.
     */
    ~MassTransport(void);



  protected:

    //! The constructor
    MassTransport(const ModelOptions& options);

    //! The initialization
    virtual void do_init(void);


    //! Parse the options from the input file
    virtual void parse_options(void);


    //! Setup the available variables
    virtual void do_setup_solution_variables(void);


    //! Solve the MassTransport equation
    virtual void do_solve(void);


    //! Print some useful information
    virtual void do_print_info(void);


    //! We need to create a physical model
    virtual PhysicalModel* create_bulk_model(const ModelOptions& options,
        const Material* mat) const;

    //! We need to create boundary condition model
    virtual PhysicalModel* create_boundary_model(const ModelOptions& options,
        const MaterialBoundary* boundary) const;


    //! We have to provide somehow our solution variables
    virtual void get_solution_secure(const Elem* elem,
        std::map<ID, std::vector<double> >& values,
        const std::vector<Point>& p);



  private:

    //! These are the known solution variables
    enum Solutions
    {
      PartialPressure,  /*!< the partial pressure */
      Concentration,    /*!< the concentration given by Henry's law */
      RelativeHumidity, /*!< the relative humidity for H2O, from the partial pressure */
      Flux,             /*!< the flux */
      Solubility,       /*!< the solubility */
      Diffusivity       /*!< the diffusion constant cm^2/s */
    };

    //! The assembly function
    void assemble(void);

    // A local helper class to be used to access assembly routine
    class MyAssembly : public TiberLinearSystem::Assembly
    {
      public:
        MyAssembly(MassTransport* obj) : _obj(obj) {};

        void assemble() override
        {
          _obj->assemble();
        }

      private:
        MassTransport *_obj;
    };

    //! The assembly object
    MyAssembly _my_assembly;


    //! The cell temperature
    double _cell_temp = 300;

    //! The molecule to be considered
    string _mass_tran_spec = "H2O"; // this can be defined in the input file

};





#endif // TC_MASSTRANSPORT_H
