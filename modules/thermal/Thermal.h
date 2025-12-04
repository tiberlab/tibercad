/*  
 * This file is part of the tiberCAD module thermal.
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
 * \file Thermal.h
 * \brief tiberCAD thermal module header.
 *
 * \note This file is part of module thermal.
 */


#ifndef _THERMAL_H_
#define _THERMAL_H_

#include "tibercad/module/SimulationInterface.h"
#include "tibercad/geom/ElementSide.h"
#include "tibercad/module/SimulationEnvironment.h"
#include "tibercad/base/tiber_dll.h"

#include "tibercad/solver/TiberLinearSystem.h"

/*!
 *
 * \brief This is an example implementation of the MyPoisson equation to
 *        help module development.
 *
 * Illustrates the basic usage of the SimulationInterface API.
 */
class TBDLLOCAL Thermal : public SimulationInterface
{

  public:

    //! Destructor
    /*!
     * We do not declare it virtual here, as we will not allow
     * to derive from this class anyway.
     */
    ~Thermal(void);

    //! We need a public static creator function
    static Thermal* create(const ModelOptions& options);


    //! The assembly function
    void assemble(void);


  protected:

    //! The initialization
    virtual void do_init(void);

    //! Parse the options from the input file
    virtual void parse_options(void){};

    //! Setup the available variables
    virtual void do_setup_solution_variables(void);

    //! Solve the MyPoisson equation
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

    //! Get a mesh independent solution variable
    virtual void get_solution_secure(std::map<ID, std::vector<double> >& values);


  private:

    // A local helper class to be used to access assembly routine
    class MyAssembly : public TiberLinearSystem::Assembly
    {
      public:
      MyAssembly(Thermal* obj) : _obj(obj) {};

      void assemble() override
      {
        _obj->assemble();
      }

      private:
      Thermal* _obj;

    };

    MyAssembly _my_assembly;

    double compute_power_dissipated();

    double compute_power_emitted();


    //! These are the known solution variables
    /*!
     * This is an enum, but we use the string representation of
     * the enum values to refer to solutions for plotting or
     * for data exchange with other modules.
     *
     * \note Do \em not use (\c INVALID_ID - 1) or the strings \c RegionIDs
     * or \c materials as they are used to plot the materials/region IDs.
     *
     * \note The name "all" is used to plot all solutions
     */
    enum Solutions
    {
      LatticeTemp,       /*!< the Lattice Temperature */
      ThermalFlux,              /*!< the thermal flux */
      HeatSource,                /*!< the HeatSource */
      ThermCond,                /*!< Thermal conductivity */
      MaxTemp                   /*!< MaxTemp */
    };

    //! The constructor
    /*!
     * Being private disables further inheritance.
     */
    Thermal(const ModelOptions& options);

    //! The maximum temperature calculated
    double _max_temperature;


};



#endif // _MYPOISSON_H_
