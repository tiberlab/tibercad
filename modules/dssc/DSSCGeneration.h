// $Id$

#include "tibercad/module/SimulationInterface.h"



class DSSCGeneration : public SimulationInterface
{


  public:

    //! Create an DSSCGeneration object
    static DSSCGeneration* create(const ModelOptions& options);


    virtual ~DSSCGeneration(void);

    
    //! Create the bulk model
    virtual PhysicalModel*
      create_bulk_model(const ModelOptions& options,
          const Material* mat) const;


    //! Create the boundary model
    virtual PhysicalModel* create_boundary_model(const ModelOptions& options,
        const Material* material_A, const Material* material_B) const;



  protected:

    //! Constructor
    DSSCGeneration(const ModelOptions& options);


    //! Initialize the equation system
    virtual void do_init(void);


    //! Solve the drift-diffusion problem.
    virtual void do_solve(void);


    //! Parse the options
    virtual void parse_options(void);


    //! Setup the available variables
    virtual void do_setup_solution_variables(void);



    /*!
     * \copydoc SimulationInterface::get_solution_secure(const Elem*,
     *  std::map<ID, std::vector<double> >&, const std::vector<Point>&)
     */
    virtual void get_solution_secure(const Elem* elem,
        std::map<ID, std::vector<double> >& values,
        const std::vector<Point>& p);



  private:

    enum Solutions
    {
      Generation,
      Distance
    };


    //! \c true if distances have been calculated
    bool _d_calculated;


    //! Direction of incident light
    libMesh::RealVectorValue _direction;


    //! Light intensity in units of 1 sun
    double _intensity;


    //! the absorption coefficient for the simple model
    double _alpha;


    //! The wavelengths
    std::vector<double> _lambda;


    //! The solar spectrum
    std::vector<double> _spectrum;


    //! The absorption spectrum
    std::vector<double> _absorption;


    //! Calculates the distances
    void _calculate_distances(void);

    //! Read the spectrum
    void _read_spectrum(void);

    //! Integrate the spectrum for a given distance
    double _integrate(double d);

};



//
// inline member functions
//




