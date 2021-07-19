// $Id$

#ifndef _TMM_H_
#define _TMM_H_

#include "SimulationInterface.h"


/*!
 * \brief Implementation of Transfer Matrix Method for electromagnetic fields
 *
 * This class implements standard Transfer Matrix Method (TMM) to calculate the
 * electromagnetic field distribution in a 1D layered structure.
 *
 * Author:
 * Contributors:
 */
class TBDLLOCAL Tmm : public SimulationInterface
{

  public:

    //! Destructor
    /*!
     * We do not declare it virtual here, as we will not allow
     * to derive from this class anyway.
     */
    ~Tmm(void);

    //! We need a public static creator function
    static Tmm* create(const ModelOptions& options);



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

    /*!
     * \brief The known solution variables
     */
    enum Solutions
    {
      EField,
      HField,
      Intensity,
      GenerationRate
    };
    
    /*!
     * \brief Constructor
     *
     * Being private disables further inheritance.
     */
    Tmm(const ModelOptions& options);

    /*!
     * \brief The wavelengths
     */
    std::vector<double> _wavelengths;

    /*!
     * \brief The incident angle
     */
    std::vector<double> _incident_angle;





};


#endif // _TMM_H_
