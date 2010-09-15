// $Id$

#ifndef _GRAYMODEL_H_
#define _GRAYMODEL_H_

#include "HeatTransportModel.h"

#include "point.h"
#include "tensor_value.h"
#include "vector_value.h"
#include "tiber_dll.h"

class Elem;


//! The base class for Poisson boundary conditions
class TBDLLOCAL GrayModel : public HeatTransportModel
{

  public:

 //! Creator function
  static GrayModel* create(const ModelOptions& options);

  //! Destructor
  ~GrayModel(void) {};


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
    GrayModel(const ModelOptions& options);

 //  struct options
//   {
    
//     double max_iter;
//     double max_error; //!< Max tollerance for self-consistent loop
//     vector<ID> spec;
//     std::vector<double> d_omega;
//     std::vector<Point> directions;
//     std::vector<Point> dir;
//     std::vector<double> theta_vec;
//     std::vector<double> phi_vec;
//     ID theta_slices;
//     ID phi_slices;
   
//   };

//  options myopts;


};



inline
PhysicalModelInterface*
GrayModel::create_new(void) const
{
  return new  GrayModel(get_options());
}

inline
GrayModel*
GrayModel::create(const ModelOptions& options)
{
  return new  GrayModel(options);
}




#endif // _GRAYMODEL_H_
