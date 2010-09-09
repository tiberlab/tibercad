// $Id$

#ifndef _GRAYMODEL_H_
#define _GRAYMODEL_H_

#include "HeatTransportModel.h"

#include "point.h"
#include "tensor_value.h"
#include "vector_value.h"
//#include "tiber_dll.h"

class Elem;


//! The base class for Poisson boundary conditions
//class TBDLLOCAL GrayModel : public HeatTransportModel
class TBDLLOCAL GrayModel : public HeatTransportModel
{

  public:

 
  //! Destructor
  ~GrayModel(void) {};

  //! Creator function
  static GrayModel* create(const ModelOptions& options);


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
  

  struct options
  {
    
    double max_iter;
    double max_error; //!< Max tollerance for self-consistent loop
    vector<ID> spec;
    std::vector<double> d_omega;
    std::vector<Point> directions;
    std::vector<Point> dir;
    std::vector<double> theta_vec;
    std::vector<double> phi_vec;
    ID theta_slices;
    ID phi_slices;
    //ID N_theta;
    //ID N_phi;

    //std::vector<int> custom_dir;
    //bool diffusive;
    //New
    //double equilibrium_energy;
    //std::string first_guess;
    //std::vector<Point> cd;
    //double s_0;
    //double t_0;
  };

 options myopts;


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
  std::cout<<"GRAY"<<std::endl;
  return new  GrayModel(options);
}




#endif // _GRAYMODEL_H_
