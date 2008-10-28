// $Id$

#ifndef _EMBRACING_H_
#define _EMBRACING_H_

#include "TypeDefs.h"

#include "point.h"
#include "linear_implicit_system.h"

#include <set>

class SimulationInterface;
class Control;
class ModelOptions;
class Elem;
class Point;



//! The embracing region of two different models
class Embracing
{

  public:

    //! The constructor
    /*!
     * \param outer the higher level simulation, typically classical
     * \param inner the lower level simulation, typically quantum
     */
    Embracing(SimulationInterface* outer, SimulationInterface* inner);


    //! The destructor
    ~Embracing(void);


    //! Initialize
    void init(const ModelOptions& options);


    //! Calculate the mixing coefficient on the embracing region
    void calculate_mixing(void);


    //! Get mixing coefficient
    /*!
     * The mixing coefficient is 0 on the outer boundary,
     * 1 on the inner boundary.
     * If \c elem is not contained in the embracing region,
     * 1.0 is returned
     */
    double get_mixing_coefficient(const Elem* elem, const Point& p);



  private:

    //! The class that handles the solve
    class LaplaceEq
      : public LinearImplicitSystem
    {

      public:

        LaplaceEq(EquationSystems& eq,
            const std::string& name, unsigned int number)
          : LinearImplicitSystem(eq, name, number) { };

        virtual void user_assembly(void) { _emb->assembly(*this); };

        void build_nodal_results(std::vector<double>& results,
            std::vector<std::string>& legend);

        Embracing* _emb;

    };


    //! The 'outer' (classical) simulation
    SimulationInterface* _outer;

    //! The 'inner' (quantum) simulation
    SimulationInterface* _inner;


    //! The embracing length
    double _lambda;


    //! \c true if the region is empty
    bool _is_empty;


    //! \c true if embracing region should be plotted
    bool _do_plot;


    //! A list of all elements that make part of this embracing region
    /*! 
     * The double corresponds to a weight indicating the distance to the boundary
     */
    std::map<const Elem*, double> _elem_list;


    //! A list of all outer boundary planes
    /*!
     * For the meaning of \em outer see Embracing::find_boundary()
     */
    std::set<ElementSide> _sides;


    //! A list of all inner boundary planes
    /*!
     * For the meaning of \em inner see Embracing::find_inner_boundary()
     */
    std::set<ElementSide> _inner_sides;


    //! The Laplace solver
    LaplaceEq* _laplace;


    //! A counter that counts the instances of this class
    static unsigned int _counter;


    //! Generates the embracing region
    void generate_embracing_region(void);


    //! Find the elements on the boundary between 'inner' and 'outer' simulation
    void find_boundary(void);

    
    //! Find the inner boundary
    /*! 
     * With \em inner boundary we intend the automatically created boundary
     * of the embracing region that does not touch the 'outer' simulation
     *
     * \note {This method assumes that the elements of the embracing region
     * are activated.}
     */
    void find_inner_boundary(void);


    //! Prepare for solving
    /*!
     * Mark all elements of the embracing region as active
     */
    void prepare_for_solve(void);


    //! Flag all elements as active
    void reactivate_all_elements(void);


    //! The matrix and rhs assembly
    void assembly(LaplaceEq& system);


    //! Plot the embracing region
    void plot(void);

};



//
// inline members
//


#endif // _EMBRACING_H_
