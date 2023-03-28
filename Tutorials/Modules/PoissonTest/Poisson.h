// $Id$

#ifndef _POISSON_H_
#define _POISSON_H_

#include "SimulationInterface.h"
#include "TiberLinearSystem.h"

#include "libmesh/enum_order.h"
#include "libmesh/enum_quadrature_type.h"

/*!
 * 
 * \brief This is an example implementation of the Poisson equation to
 *        help module development.
 *
 * Illustrates the basic usage of the SimulationInterface API.
 */
class TBDLLOCAL Poisson : public SimulationInterface
{

  public:

    //! Destructor
    /*!
     * We do not declare it virtual here, as we will not allow
     * to derive from this class anyway.
     */
    ~Poisson(void);

    //! We need a public static creator function
    static Poisson* create(const ModelOptions& options);



  protected:

    //! The initialization
    virtual void do_init(void);


    //! Parse the options from the input file
    virtual void parse_options(void);


    //! Setup the available variables
    virtual void do_setup_solution_variables(void);


    //! Solve the Poisson equation
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
      Potential,        /*!< the potential */
      ElField,            /*!< the field (negative gradient of potential) */
      Displacement,     /*!< the electric displacement */
      Polarization,
      ChargeDensity     /*!< the source (charge density) */
    };

    //! The constructor
    /*!
     * Being private disables further inheritance.
     */
    Poisson(const ModelOptions& options);

    //! The assembly function
    void assemble(void);

    // A local helper class to be used to access assembly routine
    class MyAssembly : public TiberLinearSystem::Assembly
    {
      public:
        MyAssembly(Poisson* obj) : _obj(obj) {};

        void assemble() override
        {
          _obj->assemble();
        }

      private:
        Poisson *_obj;
    };

    MyAssembly _my_assembly;


    struct Options
    {
      
      /**
       * Set the default boundary conditions
       */
      std::string default_boundary_conditions;
    
      libMeshEnums::QuadratureType quadrature_type; 

      libMeshEnums::Order integration_order;  
    };

    Options myopts;

    //! Circumcenter
    /*!
     * \param elem the element
     * \param s the side of interest (if any)
     * 
     * If s is provided, the circle will be constructed
     * such that the points of side \c s are lying on it.
     */
    libMesh::Point circumcenter(const libMesh::Elem* elem,
                                int s = -1);
    
    //! Get the geometric parameters for the coupling
    /*!
     * \param elem [in] the current element worked on
     * \param k [in] the index of the side of interest
     * \param x_i [in/out] the circumcenter of \c elem
     * \param x_k [out] the circumcenter of the neighbor
     * \param n_ik [out] the outpointing normal
     * \param x_m [out] the projection of x_i on the side
     * \param h_im [out] the distance \c x_i to \c x_m
     * \param h_mk [out] the distance \c x_m to \c x_k
     * \param A_ik [out] the volume of side \c k
     * 
     * NOTE: \c x_i can be changed in the method if \c elem
     * is not a triangle. In that case, a suitable center will
     * be calculated using three points including side \c k.
     */
    void geometry_params(const libMesh::Elem *elem,
                         unsigned int k,
                         libMesh::Point &x_i,
                         libMesh::Point &x_k,
                         libMesh::Point &n_ik,
                         libMesh::Point &x_m,
                         double &h_im, double &h_mk,
                         double &A_ik);
};





#endif // _POISSON_H_
