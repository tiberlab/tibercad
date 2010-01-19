// $Id$

#ifndef _EMBRACING_H_
#define _EMBRACING_H_

#include "ElementSide.h"

#include "point.h"
#include "linear_implicit_system.h"

#include <set>

class SimulationInterface;
class Control;
class ModelOptions;
class Elem;
class Point;



//! The embracing region of two different models
/*!
 * The default behaviour is to not calculate the mixing coefficient
 * on the embracing region. If it is needed it has to be requested
 * explicitly using the options or calling need_mixing_coeff()
 * \em before init().
 */
class Embracing
{

  public:

    //! The type of the element Map
    typedef std::map<const Elem*, double> MapType;


    //! An iterator to iterate over the elements
    class elem_iterator
    {
      public:

        //! Default constructor
        elem_iterator(void) { };

        //! Copy constructor
        elem_iterator(const elem_iterator& other)
          : _it(other._it) { };

        //! A more useful constructor
        elem_iterator(const MapType::const_iterator& it) : _it(it) { }; 

        //! Prefix increment
        elem_iterator& operator++(void)
        {
          ++_it;
          return *this;
        }

        //! Prefix increment
        elem_iterator& operator--(void)
        {
          --_it;
          return *this;
        }

        //! Assignement
        elem_iterator& operator=(const elem_iterator& rhs)
        {
          _it = rhs._it;
          return *this;
        }

        //! Comparison
        bool operator==(const elem_iterator& rhs)
        {
          return (_it == rhs._it);
        }

        //! Comparison
        bool operator!=(const elem_iterator& rhs)
        {
          return ! (*this == rhs);
        }

        //! Dereference
        const Elem* operator*(void)
        {
          return _it->first;
        }

      private:

        //! The internal iterator
        MapType::const_iterator _it;

    };


    //! The constructor
    /*!
     * \param outer the "outer" simulation,
     *  \f$\Omega_{emb}\cap\Omega_{outer} = \emptyset \f$
     * \param inner the "inner" simulation containing the
     * embracing region,
     *  \f$\Omega_{emb}\cap\Omega_{inner} = \Omega_{emb}\f$
     */
    Embracing(SimulationInterface* outer, SimulationInterface* inner);


    //! The destructor
    ~Embracing(void);


    //! Initialize
    void init(const ModelOptions& options);


    //! If we should calculate the mixing coefficients
    void need_mixing_coeff(bool need_mixing = true);


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


    //! Check if an element is in the embracing region
    bool is_in_embracing_region(const Elem* elem) const;


    //! Get an iterator to first element
    elem_iterator elem_begin(void) const;

    //! Get the past the end iterator
    elem_iterator elem_end(void) const;

    //! Find an element and return its iterator
    elem_iterator find_elem(const Elem* elem) const;


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


    //! The cutoff length
    double _cutoff;


    //! \c true if the region is empty
    bool _is_empty;


    //! \c true if embracing region should be plotted
    bool _do_plot;


    //! \c true if mixing coefficients are needed
    bool _need_mixing;


    //! A list of all elements that make part of this embracing region
    /*! 
     * The double corresponds to a weight indicating the distance to
     * the boundary.
     */
    MapType _elem_list;


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
     * With \em inner boundary we intend the automatically created
     * boundary of the embracing region that does not touch the
     * 'outer' simulation
     * \note {This method assumes that the elements of the embracing
     * region are activated.}
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

inline
Embracing::elem_iterator
Embracing::elem_begin(void) const
{
  return elem_iterator(_elem_list.begin());
}


inline
Embracing::elem_iterator
Embracing::elem_end(void) const
{
  return elem_iterator(_elem_list.end());
}


inline
Embracing::elem_iterator
Embracing::find_elem(const Elem* elem) const
{
  return elem_iterator(_elem_list.find(elem));
}


inline
bool
Embracing::is_in_embracing_region(const Elem* elem) const
{
  bool ans = false;
  if (find_elem(elem) != elem_end())
    ans = true;

  return ans;
}
 

inline
void
Embracing::need_mixing_coeff(bool need_mixing)
{
  _need_mixing = need_mixing;
}


#endif // _EMBRACING_H_
