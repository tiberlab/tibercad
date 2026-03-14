/*  
 * This file is part of the tiberCAD module pvmodule.
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
 * \file PVModule.h
 * \brief tiberCAD pvmodule module header.
 *
 * \note This file is part of module pvmodule.
 */


#ifndef TC_PVMODULE_H
#define TC_PVMODULE_H

#include "tibercad/module/SimulationInterface.h"
#include "tibercad/solver/TiberLinearSystem.h"



/*!
 * 
 * \brief A lumped-element based model for PV modules
 *
 * This model implements a Spice-based photovotaic module model.
 * It is based on a subdivision of the module into small elementary
 * cells. Each such cell is represented by an equivalent circuit,
 * and individual elementary cells are interconnected by resistors to
 * neighboring cells, representing the transport and contact layers. 
 */
class TC_DLLOCAL PVModule : public SimulationInterface
{

  public:

    //! Destructor
    ~PVModule(void);



  protected:

    //! The constructor
    explicit PVModule(const ModelOptions& options);

    //! The initialization
    virtual void do_init(void) final;


    //! Parse the options from the input file
    virtual void parse_options(void) final;


    //! Setup the available variables
    virtual void do_setup_solution_variables(void) final;


    //! Solve the PVModule equation
    virtual void do_solve(void) final;


    //! Print some useful information
    virtual void do_print_info(void) final;
	
    //! Print global data
    virtual void plot_globaldata(void) final;


    //! We need to create a physical model
    virtual PhysicalModel* create_bulk_model(const ModelOptions& options,
        const Material* mat) const final;

    //! We need to create boundary condition model
    virtual PhysicalModel* create_boundary_model(const ModelOptions& options,
        const MaterialBoundary* boundary) const final;


    //! Provide spatial solution variables
    virtual void get_solution_secure(const Elem* elem,
        std::map<ID, std::vector<double> >& values,
        const std::vector<Point>& p) final;


    //! Provide the global solutions
    virtual void get_solution_secure(std::map<ID, std::vector<double> >& values) final;

  private:

    //! These are the known solution variables
    enum Solutions
    {
      TopPotential,     /*!< the potential on top surface*/
      BottomPotential,  /*!< the potential on bottom surface */
      CellPotential,    /*!< the potential difference top - bottom*/
      CurrentDensity,   /*!< the local current density */
      ContactCurrent,   /*!< the currents at the contacts */
    };

    //! The discretization type
    enum Discretization
    {
      FEM,    /*!< FEM */
      DEC     /*!< DEC with dual grid */
    };

    //! The assembly function
    void assemble(void);

    //! The assembly function using FEM
    void assemble_fem(void);
    
    //! The assembly function using DEC on dual grid
    void assemble_dec_dual(void);

    //! Calculate current density in active nodes from the node voltages
    void calculate_current_density(void);


    // A local helper class to be used to access assembly routine
    class MyAssembly : public TiberLinearSystem::Assembly
    {
      public:
        MyAssembly(PVModule* obj) : _obj(obj) {};

        void assemble() override
        {
          _obj->assemble();
        }

      private:
        PVModule *_obj;
    };

    MyAssembly _my_assembly;

    //! The discretization to be used
    Discretization _discretization = FEM;

    std::string _spice {"ngspice"};

    //! The applied voltage
    double _voltage = 0.0;

    //! The contact current
    double _current;

    //! The node ids for ground nodes
    std::set<unsigned int> _gnd_ids;

    //! The node ids for the voltage source
    std::set<unsigned int> _src_ids;

    //! The small resistance used to connect voltage source
    double _rsource = 0.001;

    //! The small resistance used to connect ground nodes
    double _rgnd = 0.001;
	

};





#endif // TC_PVMODULE_H
