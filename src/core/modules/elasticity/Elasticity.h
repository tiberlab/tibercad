// $Id$

#ifndef _ELASTICITY_H_
#define _ELASTICITY_H_

#include "SimulationInterface.h"
#include "tensor_value.h"
#include "tensor.h"
#include "tiber_dll.h"
#include "Device.h"
/*!
 * 
 * \brief This is an example implementation of the Poisson equation to
 *        help module development.
 *
 * Illustrates the basic usage of the SimulationInterface API.
 */
class TBDLLOCAL Elasticity : public SimulationInterface
{

  public:

    //! Destructor
    /*!
     * We do not declare it virtual here, as we will not allow
     * to derive from this class anyway.
     */
    ~Elasticity(void);

    //! We need a public static creator function
    static Elasticity* create(const ModelOptions& options);



  protected:

    //! The initialization
    virtual void do_init(void);

    //! Parse the options from the input file
    virtual void parse_options(void);


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



  private:

  //!shape iteration index
  ID shape_iteration;

  //!node connection
  std::vector<unsigned short int> node_conn;

  //! A pointer to device
  Device* _device;

  //! Compute the elastic energy
  Real compute_elastic_energy(void);

  //!Get the total stress
  RealTensor get_stress(const Elem* elem, const Point& p);
  
  //!Get the internal stress
  RealTensor get_internal_stress(const Elem* elem, const Point& p);

  typedef std::map<const Elem*, RealGradient> force_vector;


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

  struct options
  {
    ID magnification;   //!< Magnification  
    ID edge;
    bool non_linear_strain;
    bool deformation;
    double shape_error; //!< Max tollerance for shape deformation
    unsigned int shape_iterations; //! Max number of shape itarations
    std::string structure_to_be_strained; //!Atomistic strucure to be strained
   };

  options myopt;


    enum Solutions
    {
      Strain,        /*!< the strain */ 
      StrainCell,        /*!< the strain */
      StrainCrystal,  /*!< the strain in the crystal system*/
      Stress,            /*!< the total stress */
      InternalStress,   /*!< the internal stress */
      InternalStrain,   /*!< the internal strain */
      StressCrystal,  /*!< the stress in the crystal system*/
      Displacement,     /*!< the displacement */
      ForceSource,      /*!< force source */
      StrainSource,     /*!< strain source */
      StressSource,     /*!< stress source */
      EnergyDensity,    /*!< elastic energy density */
      Energy             /*!< Force */
    };


  //! The total displacement
 AutoPtr<NumericVector<Number> > sol;

  force_vector internal_force;

    std::vector<unsigned int> uvar;
    //! The constructor
    /*!
     * Being private disables further inheritance.
     */
    Elasticity(const ModelOptions& options);

    //! The assembly function
    static void assemble(EquationSystems& es, const std::string& system_name);

    //! The real assembly function
    void do_assemble(EquationSystems& es, const std::string& system_name);

    //! A static pointer to this
    static Elasticity* _this;

    //! Apply the deformation
    void apply_shape_deformation();


    //! Restore the shape
    void restore_shape();


  RealTensor get_subtensor(const Tensor4DSym& C_calc,unsigned int i,unsigned  int k);
};


void
Elasticity::assemble(EquationSystems& es, const std::string& system_name)
{
  _this->do_assemble(es, system_name);
}


#endif // _ELASTICITY_H_
