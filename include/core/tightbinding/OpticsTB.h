
#ifndef _OPTICSTB_H_
#define _OPTICSTB_H_

#include "Optics.h"
#include <tensor.h>
#include "EigenvalueProblem.h"
class Mesh;
class Elem;

//! A base class of optics calculation
/*!
 * The task of this class is to calculate the spectrum at a certain given
 * k-point
 */
class OpticsTB : public Optics
{

  public:

    //! The constructor
    OpticsTB(void);

    //! The destructor
    virtual ~OpticsTB(void);

 
    static OpticsTB* create();

    void calculate_spectrum(const Mesh& Energy, double Gamma,const Tensor1& polariz, 
                          std::map<const Elem*, double>& spectrum);

  protected:

    virtual void do_init(void);

    virtual void do_solve(void);

    virtual void do_plot(void);
    
    virtual void parse_options (void);

  private:

    //! Mesh for spectrum [eV];  (should go in Optics.h)
    Mesh* _energy_mesh;

    //!  momentum matrix elements (should go in Optics.h)
    /*!
      Px_matrix[j][k]:  j - initial state; k - final state
    */    
    std::vector <std::vector <Complex> >    Px_matrix;
    std::vector <std::vector <Complex> >    Py_matrix;
    std::vector <std::vector <Complex> >    Pz_matrix;

    //!numbers of eigensates that are considered as intial states for optical transition
    std::vector<unsigned int> _initial_eigen_state_numbers;
    
    
    //!numbers of eigensates that are considered as final states for optical transition
    std::vector<unsigned int> _final_eigen_state_numbers;
    

    //!type of particle for the intial states 
    std::string _initial_state_particle;


    //!type of particle for the final states 
    std::string _final_state_particle;
    


    //!pointer to the eigenvalue solver for initial states
    EigenvalueProblem* _initial_state_model;
    
    
    //!pointer to the eigenvalue solver for final states
    EigenvalueProblem* _final_state_model;

    //! initial states
    std::vector<EigenvalueProblem::eigen_problem_solution> _i_states;    

    //! final states
    std::vector<EigenvalueProblem::eigen_problem_solution> _f_states;       

    void calculate_P_matrix_elements( );

};


inline OpticsTB* OpticsTB::create()
{
  return (new OpticsTB);
}


#endif // _OPTICSTB_H_
