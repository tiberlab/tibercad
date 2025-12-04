/*  
 * This file is part of the tiberCAD module wateringress.
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
 * \file WaterIngress.h
 * \brief tiberCAD wateringress module header.
 *
 * \note This file is part of module wateringress.
 */


#ifndef _WATERINGRESS_H_
#define _WATERINGRESS_H_

#include "tibercad/module/SimulationInterface.h"
#include "tibercad/solver/TiberLinearSystem.h"


/*!
 * 
 * \brief This is a simple implementation of the Fick's laws for
 *        mmodeling water ingress e.g. into solar cells.
 *
 * The implementation uses the partial pressure as primary
 * variable, and assumes the solubility to be piecewise constant.
 */
class TBDLLOCAL WaterIngress : public SimulationInterface
{

  public:

    //! Destructor
    /*!
     * We do not declare it virtual here, as we will not allow
     * to derive from this class anyway.
     */
    ~WaterIngress(void);

    //! We need a public static creator function
    static WaterIngress* create(const ModelOptions& options);



  protected:

    //! The initialization
    virtual void do_init(void);


    //! Parse the options from the input file
    virtual void parse_options(void);


    //! Setup the available variables
    virtual void do_setup_solution_variables(void);


    //! Solve the WaterIngress equation
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
      PartialPressure,  /*!< the partial water pressure */
      Concentration,    /*!< the water concentration given by Henry's law */
      RelativeHumidity, /*!< the relative humidity, from the partial pressure */
      Flux,             /*!< the flux */
      Solubility,       /*!< the solubility */
      Diffusivity       /*!< the diffusion constant cm^2/s */
    };

    //! The constructor
    /*!
     * Being private disables further inheritance.
     */
    WaterIngress(const ModelOptions& options);

    //! The assembly function
    void assemble(void);

    // A local helper class to be used to access assembly routine
    class MyAssembly : public TiberLinearSystem::Assembly
    {
      public:
        MyAssembly(WaterIngress* obj) : _obj(obj) {};

        void assemble() override
        {
          _obj->assemble();
        }

      private:
        WaterIngress *_obj;
    };

    //! The assembly object
    MyAssembly _my_assembly;


    //! The cell temperature
    double _cell_temp = 300;

};





#endif // _WATERINGRESS_H_
