
#ifndef _OPTICSTB_H_
#define _OPTICSTB_H_

#include "Optics.h"
#include <tensor.h>
#include "EigenvalueProblem.h"

//! A base class of optics calculation
/*!
 * The task of this class is to calculate the spectrum at a certain given
 * k-point
 */
class OpticsTB : public Optics
{

  public:

    //! The constructor
    OpticsTB(const ModelOptions& options);

    //! The destructor
    virtual ~OpticsTB(void);

 
    static OpticsTB* create(const ModelOptions& options);


  protected:

    virtual void do_init(void);

    virtual void do_assemble(const ModelOptions& options);

    //! Assemble the P-matrix and compute its matrix elements.
    virtual void do_compute_matrix_elements(void);

    
    virtual void calculate_matrix_bulk(void){};


  private:


    //! checks that states for optics are really there
    void check_states(void);

};


inline OpticsTB* OpticsTB::create(const ModelOptions& options)
{
  return (new OpticsTB(options));
}


#endif // _OPTICSTB_H_
