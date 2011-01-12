// $Id$

#include "SimulationInterface.h"


class DSSCGeneration : public SimulationInterface
{

  public:

    enum Solutions
    {
      Generation
    };

    //! Create an DSSCGeneration object
    static DSSCGeneration* create(const ModelOptions& options);


    virtual ~DSSCGeneration(void) {};

    
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

    //! \c true if distances have been calculated
    bool _d_calculated;


    //! Direction of incident light
    RealVectorValue _direction;


    //! Light intensity
    double _intensity;


    //! Calculates the distances
    void _calculate_distances(void);

};



//
// inline member functions
//

inline
DSSCGeneration*
DSSCGeneration::create(const ModelOptions& options)
{
  return new DSSCGeneration(options);
}


