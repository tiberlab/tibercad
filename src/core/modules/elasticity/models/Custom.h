// $Id: Clamp.h 2124 2010-10-22 14:00:17Z gromano $

#ifndef _Custom_H_
#define _Custom_H_

#include "ElasticityBoundaryModel.h"

#include "point.h"
#include "tensor_value.h"
#include "vector_value.h"



class Elem;


//! The base class for Poisson boundary conditions
class Custom : public ElasticityBoundaryModel
{

  public:

    //! Destructor
    ~Custom(void) {};

    //! Creator function
    static Custom* create(const ModelOptions& options);


    //! Calculate for a point on the given side
    virtual void calculate(const Elem* elem, unsigned int side,
			   const Point& point){};


  protected:

    //! Initialize
  virtual void do_init(void);

    /* In some cases it might be useful to reimplement this: */
    // virtual void do_init_interface(const PhysicalModelInterface* comp_A,
    //         const PhysicalModelInterface* comp_B);


    /* This is not used here: */
    // virtual void read_database(void);


    /* We do not use this here: */
    // virtual void read_interface_database(void);


    //! Create a new object of the same type
    virtual PhysicalModelInterface* create_new(void) const;


  private:

    //! Constructor
    Custom(const ModelOptions& options);

   
};



inline
Custom::Custom(const ModelOptions& options) :
  ElasticityBoundaryModel(options)
{
}



inline
Custom*
Custom::create(const ModelOptions& options)
{
  return new Custom(options);
}



inline
PhysicalModelInterface*
Custom::create_new(void) const
{
  return new Custom(get_options());
}

#endif // _POISSONDIRICHLET_H_
