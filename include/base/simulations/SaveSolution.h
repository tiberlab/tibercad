// $Id: SaveSolution.h 3735 2014-01-23 09:04:13Z maufder $

#ifndef _SAVESOLUTION_H_
#define _SAVESOLUTION_H_

#include "SimulationInterface.h"

#include <set>
#include <vector>
#include <string>



//! Make a sweep of a sweepable variable
/*!
 * This class takes a variable (cf. Variable class)
 * and performs a sweep with it (e.g. contact voltage). For each
 * sweep step it solves the given simulation(s) according to
 * the following flow chart:
 * \image html Sweep_flowchart.jpg
 * \image latex Sweep_flowchart.eps
 *    "Flow chart for a parameter sweep" width=10cm
 */
class TBDLLOCAL SaveSolution : public SimulationInterface
{

  public:
    
    //! Destructor
    virtual ~SaveSolution(void);

    //! Create a Sweep object
    static SaveSolution* create(const ModelOptions& options);


    //! Get the list of inner simulations
    const std::vector<SimulationInterface*>&
      get_simulations(void) const;


  protected:

    //! Constructor
    SaveSolution(const ModelOptions& options);

    /*! \copydoc SimulationInterface::do_init() */
    virtual void do_init(void);


    /*! \copydoc SimulationInterface::do_solve()
     *
     * solve() of the associated SimulationInterface objects will
     * be called repeatedly for all sweep values
     */
    virtual void do_solve(void);

    //! Abused here for the inheritance of variables
    virtual void do_print_info(void);

    //! Get global quantities from all sub-simulations
    virtual void get_solution_secure(const libMesh::Elem* elem, 
                                           std::map<ID, std::vector<double> >& values,
                                     const std::vector<libMesh::Point>& points);

    virtual void parse_options(void);

    virtual void do_setup_solution_variables(void);


  private:

    //! Build the equation system
    void build_equation_system(void);

    //! The simulations for which we save the solution
    SimulationInterface* _simulation;

    //! List of solution IDs
    std::vector<ID> _solution_ids;

    //! A map storing the saved solutions
    std::map<ID, std::vector<double>> solution;


};


//
// inline methods
//


inline
SaveSolution::SaveSolution(const ModelOptions& options)
  : SimulationInterface(options)
{

}


inline
SaveSolution*
SaveSolution::create(const ModelOptions& options)
{
  return new SaveSolution(options);
}


#endif // _SWEEP_H_
