// $Id: ThermalBalance.h 2457 2011-03-06 23:52:12Z gromano $

#ifndef _THERMAL_H_
#define _THERMAL_H_

#include "SimulationInterface.h"
#include "ElementSide.h"
#include "SimulationEnvironment.h"
#include "tiber_dll.h"

class TiberLinearSystem;

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

   /*! \copydoc SimulationInterface::do_get_solution_vector() */
    virtual NumericVector<double>& do_get_solution_vector(void);

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

  //!node connection
  std::vector<unsigned short int> node_conn;

  double compute_power_dissipated();

  double compute_power_emitted();

  //TiberLinearSystem*  system;

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

    //! The assembly function
    static void assemble(EquationSystems& es, const std::string& system_name);

    //! The real assembly function
    void do_assemble(EquationSystems& es, const std::string& system_name);

    //! A static pointer to this
    static Thermal* _this;



};


inline
void
Thermal::assemble(EquationSystems& es, const std::string& system_name)
{
  _this->do_assemble(es, system_name);
}


#endif // _MYPOISSON_H_
