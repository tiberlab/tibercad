// $Id: PVModule.h 4391 2017-04-07 11:16:58Z pecchia $

#ifndef _PVMODULE_H_
#define _PVMODULE_H_

#include "SimulationInterface.h"
#include "TiberLinearSystem.h"


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
class TBDLLOCAL PVModule : public SimulationInterface
{

  public:

    //! Destructor
    ~PVModule(void);

    //! We need a public static creator function
    static PVModule* create(const ModelOptions& options);



  protected:

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


    //! We need to create a physical model
    virtual PhysicalModel* create_bulk_model(const ModelOptions& options,
        const Material* mat) const final;

    //! We need to create boundary condition model
    virtual PhysicalModel* create_boundary_model(const ModelOptions& options,
        const MaterialBoundary* boundary) const final;


    //! We have to provide somehow our solution variables
    virtual void get_solution_secure(const Elem* elem,
        std::map<ID, std::vector<double> >& values,
        const std::vector<Point>& p) final;



  private:

    //! These are the known solution variables
    enum Solutions
    {
      TopPotential,     /*!< the potential on top surface*/
      BottomPotential,  /*!< the potential on bottom surface */
      CurrentDensity,   /*!< the local current density */
      ContactCurrent,   /*!< the currents at the contacts */
    };

    //! The constructor
    /*!
     * Being private disables further inheritance.
     */
    PVModule(const ModelOptions& options);

    //! The assembly function
    void assemble(void);

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

    std::string _spice {"ngspice"};
	
	std::vector<double> _jv_ref_v;
	std::vector<double> _jv_ref_j;

};





#endif // _PVMODULE_H_
