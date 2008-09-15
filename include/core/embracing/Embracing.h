// $Id$

#ifndef _EMBRACING_H_
#define _EMBRACING_H_

#include "TypeDefs.h"

#include "point.h"

#include <set>

class SimulationInterface;
class ModelOptions;
class Elem;
class Pnt;



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
     * some explanation...
     */
    double get_mixing_coefficient(const Elem*, const Point& p);


    void plot(void);



  private:

    //! The 'outer' (classical) simulation
    SimulationInterface* _outer;

    //! The 'inner' (quantum) simulation
    SimulationInterface* _inner;


    //! The embracing length
    double _lambda;


    //! \c true if the region is empty
    bool _is_empty;


    //! A list of all elements that make part of this embracing region
    /*! 
     * The double corresponds to a weight indicating the distance to the boundary
     */
    std::map<const Elem*, double> _elem_list;


    //! A list of all Planes
    std::set<ElementSide> _sides;


    //! Generates the embracing region
    void generate_embracing_region(void);


    //! Find the boundary elements
    void find_boundary(void);


    //! Prepare for solving
    /*!
     * Mark all elements of the embracing region as active
     */
    void prepare_for_solve(void);


    //! Flag all elements as active
    void reactivate_all_elements(void);
};





#endif // _EMBRACING_H_
