/*  
 * This file is part of the tiberCAD module degradation.
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
 * \file Degradation.h
 * \brief tiberCAD degradation module header.
 *
 * \note This file is part of module degradation.
 */


#ifndef TC_DEGRADATION_H
#define TC_DEGRADATION_H

#include "tibercad/module/SimulationInterface.h"
#include "tibercad/module/SolutionProvider.h"    
#include "tibercad/solver/TiberLinearSystem.h"

/*!
 * 
 * \brief This is a simple implementation of the kinetic model for the
 *        degradation of perovskite absorbers into solar cells.
 *
 * The implementation uses the wateringress model (partial pressure as primary
 * variable, and assumes the solubility to be piecewise constant).
 */
class TC_DLLOCAL Degradation : public SimulationInterface
{

  public:

    //! Destructor
    /*!
     * We do not declare it virtual here, as we will not allow
     * to derive from this class anyway.
     */
    ~Degradation(void);



  protected:

    //! The constructor
    Degradation(const ModelOptions& options);

    //! The initialization
    virtual void do_init(void);


    //! Parse the options from the input file
    virtual void parse_options(void);


    //! Setup the available variables
    virtual void do_setup_solution_variables(void);


    //! Solve the Degradation equation
    virtual void do_solve(void);


    virtual double do_eval_deg_rate(const libMesh::Elem* elem,
        const libMesh::Point& p);


    //! Print some useful information
    virtual void do_print_info(void);


    //! We have to provide somehow our solution variables
    virtual void get_solution_secure(const Elem* elem,
        std::map<ID, std::vector<double> >& values,
        const std::vector<Point>& p);



  private:

    //! These are the known solution variables
    enum Solutions
    {
      PerovskiteFraction  /*!< the fraction of active (i.e., absorbing) perovskite */
    };

    //! The assembly function
    void assemble(void);

    // A local helper class to be used to access assembly routine
    class MyAssembly : public TiberLinearSystem::Assembly
    {
      public:
        MyAssembly(Degradation* obj) : _obj(obj) {};

        void assemble() override
        {
          _obj->assemble();
        }

      private:
        Degradation *_obj;
    };

    //! The assembly object
    MyAssembly _my_assembly;


    //! The cell temperature
    double _cell_temp = 300; // [K]
    
    //! The photon flux generating carriers
    double _el_den = 1.5e21; // [(ph. m^-2 s^-1)^0.72]^0.72 when computing reaction rate 

    //! From where to get humidity
    SolutionProvider _humidity_model;

    //! From where to get oxygen
    SolutionProvider _oxygen_model;

};





#endif // TC_DEGRADATION_H



