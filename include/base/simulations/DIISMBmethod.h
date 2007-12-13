// $Id$
#ifndef _DIISMBMETHOD_H_
#define _DIISMBMETHOD_H_

#include "ModifiedBroyden.h"

//! Based on the paper Journal of Chemical Physics v.108 p.4426 (Appendix B)
class DIISMBmethod: public ModifiedBroyden
{
 public:

  //! The constructor
  DIISMBmethod(void) {};


  //! The destructor
  virtual ~DIISMBmethod(void) {};

  //! Create a Sweep object
  static DIISMBmethod* create(void);
  

 protected:

  //! DIISMB step: \f$ X_{n+1} = X_{n} + \alpha (s_{nn} - 1) F(X_n) + \alpha \sum_{i = 2}^{n-1} s_{ni}F(X_i)   \f$
  virtual void calculate_new_solution();
  
  //! calculate \f$ \tilde a \f$ matrix (see below)
  /*!
    \f$
    \tilde a_{n-1,i} = \mu_{n-1,n-1} (\lambda_{ni} - \lambda_{n-1,i} ) (i = 1,\ldots, n-2),
    \tilde a_{n-1,n-1} = 1 
    \f$
  */ 
  virtual void calculate_tilde_a(void);

 private:
  
  //! calcualte \f$ \lambda \f$ matrix (see below)
  /*!
    \f$
    \lambda_{i,n-1} = \mu_{n-1,n-1}(\kappa_{n,i} - \kappa_{n-1,i}) (i = 1,\ldots, n-1),\\
    \lambda_{n,i} = \mu_{i,i} (\kappa_{n,i+1}- \kappa_{n,i}) (i = 1,\ldots, n - 1)
    \f$
  */
  void calculate_lambda_matrix(void);

  //! \f$ \mbox{\boldmath$\eta$} = \mbox{\boldmath$\lambda$} \cdot \mbox{\boldmath$\omega$} \cdot \mbox{\boldmath$\zeta$} \f$
  void calculate_eta_matrix(void);


  //! calculate \f$ q \f$ matrix (see below)
  /*!
    \f$
    q_{n-1,i} = \sum_{j=2}^{n-1} s_{n-1,j}\kappa_{jn} (i = 1,\ldots,n-1),\\
    q_{i,n} = \sum_{j=2}^{i} s_{i,j}\kappa_{jn} (i = 1,\ldots,n-1)
    \f$
  */
  void calculate_q_matrix(void);

  //! calculate \f$ \sigma \f$ matrix (see below)
  /*!
    \f$
    \sigma_{n-1,i} = \alpha(q_{n-1,i} - \kappa_{n,i}) (i = 1,\ldots, n-1),\\
    \sigma_{i,n}   = \alpha(q_{i,n} - \kappa_{i+1,n}) (i = 1,\ldots, n-1)
    \f$
   */
  void calculate_sigma_matrix(void);

  //!calculate \f$ \nu \f$ matrix (see below)
  /*!
    \f$
    \nu_{i,n-1} = \nu_{n-1,i} = \alpha^2 \left(  
    \kappa_{n, i+1} - q_{n-1, i+1} - q_{i,n} + \sum_{j=2}^i s_{i,j}q_{n-1,j}
\right)  (i = 1, \ldots, n-1)
    \f$
  */
  void calculate_nu_matrix(void);

  
  //! \f$ \epsilon_{(m)} = \alpha^2 \kappa_{(m)} +  \alpha \left( \eta_{(m)} \cdot \sigma_{m} +  (\eta_{(m)} \cdot \sigma_{m})^T \right) + \eta_{(m)} \cdot \nu * \eta_{(m)^T}.\f$
  void calculate_DIIS_matrix(void);


  
  void calculate_c_m(void);

  void calculate_o_m_and_phi(void);

  void calculate_s_matrix(void);


  //! algorithm of the  DIISMB method
  virtual void do_step(void);


  NEWMAT::Matrix    _s, _lambda, _eta, _q, _sigma, _o_m, _eta_m;

  NEWMAT::SymmetricMatrix _nu; 

  NEWMAT::Matrix _eps;

  NEWMAT::ColumnVector  _c_m, _phi;

  

};

inline
DIISMBmethod* DIISMBmethod::create()
{
  return new DIISMBmethod();
}




#endif
