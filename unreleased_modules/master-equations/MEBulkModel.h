
#ifndef _MEBULKMODEL_H_
#define _MEBULKMODEL_H_



#include "MasterEquationsProperties.h"

//#include "TemperatureInterface.h"
//#include "StrainInterface.h"
#include "SimulationOptions.h"
#include "MasterEquationsDefs.h"
#include "TiberCad.h"
#include "Constants.h"
#include "TypeDefs.h"

#include "tensor.h"

#include "tensor_value.h"
#include "vector_value.h"
#include "point.h"


#include <vector>
#include <set>
#include <map>


// forward declarations
class Elem;
//class Dopant;
//class Trap;
class SimulationInterface;
class MasterEquations;
//class MobilityModelInterface;
//class ThermoelectricPower;
//class PolarizationModel;


//! The base class for all drift-diffusion related semiconductor models
/*!
 * \note { it is not yet possible to have twice the same recombination model.
 * Trying to add to identical models will result in a memory leak. This will
 * be corrected in future. }
 */
class MEBulkModel : public MasterEquationsProperties
{

  public:

    //! A default (empty) destructor.
    virtual ~MEBulkModel(void);


    //! Create a named Master Equations Bulk model
    /*!
     * The model is created according to the given model name.
     * If it is not known, the NULL pointer is returned.
     *
     * \param name the model name
     * \param options the options as given in the input file
     * \return a pointer to the newly created object
     */
    static MEBulkModel* create(const std::string& name,
        const Material* mat,
        const ModelOptions& options = ModelOptions());



    //! (Re-)Initialize for the given element
    /*!
     * \c reinit() calls \c prepare_element_data() which needs to be
     * implemented in derived classes
     */
    virtual void do_reinit(const Elem* elem);



    //! Calculates the equilibrium properties.
    /*!
     *
     * This method has to be called before calculating anything
     * Call this method from a derived one.
     *
     * \pre { \c reinit() has to be called before }
     * \post { all equilibrium properties are accessible without
     *  explicitly calling \c calculate_all() }
     */
    virtual void calculate_equilibrium_properties(void);



    //! Get the all nodal temperatures for a given element
    std::vector<double>& get_temperature_at_nodes(void);




  protected:



    //! The empty constructor.
    MEBulkModel(const ModelOptions& options);


    /*! \copydoc PhysicalModel::read_database() */
    virtual void read_database(void);


    //! Initialize this model
    /*!
     * This reads the database and calls init for all submodels
     * A derived class which reimplements this method has to call
     * explicitly the one of this class!
     */
    virtual void do_init(void);


    //! This method gets called from reinit()
    /*!
     * It can be used to setup data that is constant in an element, e.g.
     * strain related stuff, band edges.
     * This method can be used overiden by derived classes.
     */
    virtual void prepare_element_data(void) {};


    //! \copydoc PhysicalModel::do_print_info(void)
    virtual void do_print_info(void);



    //! Get the temperature interface
    TemperatureInterface& get_temperature_interface(void);


    //! Tells if we are doing equilibrium calculation
    bool has_solution(void) const;


    //! Set the intrinsic Fermi level
    /*!
     * Calculates the associated equilibrium densities
     */
    void set_equilibrium_properties(double Ef);


  private:


    //! The copy constructor is disabled
    MEBulkModel(const MEBulkModel& rhs);


    //! The assignment operator is disabled
    MEBulkModel& operator=(const MEBulkModel& rhs);


    //! Parse the model options
    void parse_options(void);


    //! The nodal lattice temperature
    std::vector<double> _nodal_lattice_vt;


};





#endif /* _MEBULKMODEL_H_ */
