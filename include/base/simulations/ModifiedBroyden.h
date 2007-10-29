
#ifndef _MBSELFCONSISTENT_H_
#define _MBSELFCONSISTENT_H_
#include "newmat.h"
#include "newmatap.h"
#include "SelfconsistentSolver.h"


class ModifiedBroyden:  public SelfconsistentSolver
{

 public:

  //! The constructor
  ModifiedBroyden(void);
  
  //! The destructor
  virtual ~ModifiedBroyden(void);

  //! Create a Sweep object
  static ModifiedBroyden* create(void);


 protected:

  /*! \copydoc SimulationInterface::do_init() */
  virtual void do_init(void);


  /*! \copydoc SimulationInterface::do_solve() */
  virtual void do_solve(void);

  
  /*! \copydoc SimulationInterface::parse_options() */
  virtual void parse_options(void);
  
  /*! \copydoc SimulationInterface::build_integrated_quantities() */
  virtual void build_integrated_quantities(
					   const std::set<std::string>& names,
					   std::vector<double>& values);

  
  /*! \copydoc SimulationInterface::build_integrated_quantities() */
  virtual void build_integrated_quantities_description(
						       const std::set<std::string>& names,
						       std::vector<std::string>& legend,
						       std::vector<std::string>& description);

  
 private:

  
  void init_X(void);


  void do_iteration(void);


  //!gets the calculated quantity
  void get_X_from_iteration(void);

  //!passes X to the iteration
  void pass_X_to_iteration(void);

  /*! \f$ F(x) = x - x' \f$  */
  void evaluate_F( );


  /*! \f$ \kappa_{in} = \kappa_{ni} =  F(x_n)^T \cdot F(x_i) \f$ */
  void calculate_kappa_matrix(void); 


  /*! \f$ x_{n+1} = x_{n} + \alpha * (p_{nn} - 1)F(x_n) + \alpha \sum_{i=2}^{n-1} p_{ni}F(X_i) \f$ */
  void calculate_new_solution();

  /*! \f$ 

    \mu_{i,n-1} = \mu_{n-1,i} = 0 (i = 1,...,n-2), \\
    \mu_{n-1, n-1} = (\kappa_{nn} - 2\kappa_{n,n-1} + k_{n-1, n -1 })^{-1/2}; \\

    w_{ij} = \omega_{i}\delta_{ij}.

    \f$ 
  */
  void calculate_mu_and_w_diag(void);

  /*!
    \f$
    \tilde a_{n-1,i} = \lambda_i (i = 1,...,n-2)
    \tilde a_{n-1,n-1} = 1
    \f$
  */
  void calculate_tilde_a(void);

  /*!
  \f$ \lambda_i = \mu_{ii} (\kappa_{n,i+1} -\kappa_{ni} ) (i = 1,...,n-1) \f$
  
  */
  void calculate_lambda(void);

  /*! 
    \f$ a_{ij} = a_{ji} = w_{ii} w_{jj} \tilde a_{ij} (i \ge j = 1,...,n-1) \f$
  */
  void calculate_a(void);


  /*!
    \f$ {\mathbf \beta} = (\omega_0^2 I + {\bf a})^{-1}  \f$
   */
  void calcvulate_beta_matrix(void);

  /*!
    \f$ {\mathbf \zeta}^{(n)} = {\bf \beta} \cdot {\bf w} \cdot {\bf \mu} + \omega_0^2 {\bf \beta'} \cdot {\bf \zeta}^{(n-1)} \f$
   */
  void calculate_zeta_matrix(void);

  /*!
    \f$ {\bf \eta } = {\bf w} \cdot {\bf \zeta}^{(n)} \cdot {\bf \lambda} \f$
  */
  void calculate_eta_vector(void);

  /*!
    \f$ p_{ni} = \eta_{ni} - \sum_{j=1}^{n-1} \eta_j p_{ji} (i = 2,...,n)\f$
  */
  void calculate_p_matrix(void);
    

  double evaluate_difference();


  NEWMAT::Matrix  _kappa, _beta, _p, _zeta, _a, _a_tilde;

  NEWMAT::ColumnVector _mu, _w, _lambda;


  //! \f$ F_i = X - X'  \f$
  std::vector < std::vector <double> > _F;


  //! current solution old
  std::vector <double> _X;
 
  //! current solution old
  std::vector <double> _X1;


  //! current iteration number
  unsigned int  _it_number;



  //! size of _X
  unsigned int _vector_size;

}; 


inline ModifiedBroyden* ModifiedBroyden::create(void)
{
  return new ModifiedBroyden();
}


#endif //_MBSELFCONSISTENT_H_
