#ifndef _AUTOMATICKEATING_H_
#define _AUTOMATICKEATING_H_

#include "Keating.h"

namespace libMesh
{
  template<typename T> class DenseMatrix;
  template<typename T> class DenseVector;
}

//! User defined Keating model parameters
class TBDLLOCAL AutomaticKeating : public Keating
{
public:

  //! Destructor
  virtual ~AutomaticKeating(void) {};

  //! Creator function
  static AutomaticKeating* create(const ModelOptions& options);

protected:

  //! Assign value to parameters
  virtual void do_init(void);

  //! Print some info
  virtual void do_print_info(void);

private:

  AutomaticKeating(const ModelOptions& options);

  //! Calculate zb alpha from stiffness constants
  void calculate_zb_alpha(void);

  //! Calculate zb beta from stiffness constants
  void calculate_zb_beta(void);

  //! Optimize wz alpha and beta to fit bulk elastic constants
  void calculate_wz_params(void);

  //! Calculate elastic constants residual from Keating parameters
  /*!
   * Ordering is C11, C33, C12, C13, C44, C66
   */
  void residual_wz(double a, double c, double u,
      const libMesh::DenseVector<double>& keating,
      libMesh::DenseVector<double>& residual,
      libMesh::DenseMatrix<double>& gradients);

  //! Get needed information from material database
  /*
   * If alpha and beta are found here, then they will be used.
   * Otherwise they will be calculated from the stiffness constant.
   * Note that the first method is not safe for alloys, as the stiffness
   * constant follow Vegard law while alpha and beta don't.
   */
  void parse_zb_database(void);

  //! Get needed information from material database
  /*
   * If alpha and beta are found here, then they will be used.
   * Calculation from stiffness constants is still not supported for wurtzite.
   */
  void parse_wz_database(void);

  double _c11;
  double _c12;
  double _c44;
  double _c13;
  double _c33;

  //! Weights for the error functional
  /*!
   * Order is: C11, C12, C33, C13, C44, C66
   * Default:  1,   1,   2,   2,   0.5, 0
   */
  std::vector<double> _weights;


  //! Whether we should use four parameters for wurtzite
  bool _use_four_wz_params;

};

inline
AutomaticKeating*
AutomaticKeating::create(const ModelOptions& options)
{

     return new AutomaticKeating(options);

}



#endif
