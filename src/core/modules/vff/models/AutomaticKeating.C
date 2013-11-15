#include "AutomaticKeating.h"
#include "Database.h"
#include "Messages.h"
#include "RuntimeException.h"

#include "TiberModule.h"

#include "dense_matrix.h"
#include "dense_vector.h"

#include "fadiff.h"

using namespace std;
using namespace fadbad;


// an anonymous namespace
namespace
{

  template<typename valueT>
  inline void keating_wz(double a, double c, double u,
      const vector<valueT>& x, vector<valueT>& y)
  {

    double sq3 = std::sqrt(3.0);

    double v = 1 - 2*u;
    double w = 1 - 4*u;
    double r0 = std::sqrt(3 * c*c * v*v + 4 * a*a) / (2 * sq3);
    double r01 = c * u;

    valueT A = x[0] / (r0*r0);
    valueT B = x[1] / (r0*r0);
    int shift = (x.size() == 4) ? 2 : 0;
    valueT A1 = x[shift] / (r01*r01);
    valueT B1 = x[shift + 1] / (r01*r01);

    double u2 = u*u;
    double v2 = v*v;
    double w2 = w*w;


    y.resize(6);

    // C11
    valueT D1 = (2*A + B) * ( 8*A1*u2 + 6*(A+2*B)*v2 + 3*B1*w2 );
    y[0] = a*a / (2*sq3*c) *
      ((4*A*A + 13*A*B + B*B) * (8*A1*u2 + 3*B1*w2) + 162*A*B*(A + B)*v2) / D1;

    // C12
    y[1] = a*a / (2*sq3*c) *
      (A-B)*((4*A-B)*(8*A1*u2 + 3*B1*w2) + 54*A*B*v2 ) / D1;

    // C33
    valueT D2 = 8*A1*u2 + 6*(A+2*B)*v2 + 3*B1*w2;
    y[2] = 3*sq3*c*c*c / (4*a*a) *
      ((A+2*B) * (8*A1*u2*v2 + 3*B1*v2*v2) + 16*A1*B1*u2*u2  ) / D2;

    // C13
    y[3] = sq3*c / 2  *
      (A-B)*(8*A1*u2 + 3*B1*v*w)*v / D2;


    // C44
    y[4] = sq3*c / 4 *
        ((2*A + B) * B1 ) / (2*A + B + B1);

    // C66
    y[5] = 0.5 * (y[0] - y[1]);

  }
}



AutomaticKeating::AutomaticKeating(const ModelOptions& options):
  Keating(options),
  _c11(0.0),
  _c12(0.0),
  _c44(0.0),
  _c13(0.0),
  _c33(0.0),
  _weights({1, 1, 2, 2, 0.5, 0}),
  _use_four_wz_params(false)
{

}

void
AutomaticKeating::do_init(void)
{
  Keating::do_init();

  std::string warning1("Calculating keating parameters for alloy with Vegard's law. This is not safe.");
  if (get_material()->get_structure() == "zb")
  {
    parse_zb_database();
    if ((alpha_0() == 0.0) || (alpha_1() == 0.0))
    {
      if (get_material()->is_alloy())
        Messages::warning(warning1);
      calculate_zb_alpha();
    }
    if ((beta_0() == 0.0) || (beta_1() == 0.0))
    {
      if (get_material()->is_alloy())
        Messages::warning(warning1);
      calculate_zb_beta();
    }
  }

  if (get_material()->get_structure() == "wz")
  {
    parse_wz_database();
    get_option("weights", _weights);
    _use_four_wz_params = get_option("use_four_wz_parameters", _use_four_wz_params);
    if ((alpha_0() == 0.0) || (alpha_1() == 0.0) ||
        (beta_0() == 0.0) || (beta_1() == 0.0))
    {
      calculate_wz_params();
    }
  }

}

void
AutomaticKeating::parse_zb_database(void)
{
  double alpha, beta;
  const Database& db = get_database();
  db.set_section("keating");

  alpha = db.get("alpha", 0.0, false);
  alpha_0() = alpha;
  alpha_1() = alpha;
  beta = db.get("beta", 0.0, false);
  beta_0() = beta;
  beta_1() = beta;

  db.set_section("elasticity");
  _c11 = db.get("C11", 0.0, true);
  _c12 = db.get("C12", 0.0, true);
  _c44 = db.get("C44", 0.0, true);


}


void
AutomaticKeating::parse_wz_database(void)
{
  double alpha, beta;
  const Database& db = get_database();
  db.set_section("keating");

  alpha = db.get("alpha", 0.0, false);
  alpha_0() = alpha;
  alpha_1() = alpha;
  beta = db.get("beta", 0.0, false);
  beta_0() = beta;
  beta_1() = beta;
  alpha = db.get("alpha_0", 0.0, false);
  alpha_0() = alpha;
  alpha = db.get("alpha_1", 0.0, false);
  alpha_1() = alpha;
  beta = db.get("beta_0", 0.0, false);
  beta_0() = beta;
  beta = db.get("beta_1", 0.0, false);
  beta_1() = beta;

  db.set_section("elasticity");
  _c11 = db.get("C11", 0.0, true);
  _c12 = db.get("C12", 0.0, true);
  _c44 = db.get("C44", 0.0, true);
  _c13 = db.get("C13", 0.0, true);
  _c33 = db.get("C33", 0.0, true);


}


void
AutomaticKeating::calculate_wz_params(void)
{
  double a = get_a();
  double c = get_c();
  double u = get_u();

  const int nvar = _use_four_wz_params ? 4 : 2;

  DenseVector<double> keating(nvar);
  keating(0) = keating(2 % nvar) = 80;
  keating(1) = keating(3 % nvar) = 20;


  // the database elastic constants
  //DenseVector<double> c_ref;
  //c_ref.get_values() = {_c11, _c12, _c33, _c13, _c44, 0.5 * (_c11 - _c12)};

  // for the calculated elastic constants minus "exact" values
  DenseVector<double> residual;

  // size will be provided
  DenseMatrix<double> grad;

  // We use a least squares approach to find the alpha and beta
  // for this we write f = 1/2*sum_i(r_i*r_i)
  // with r_i = c\_calc_i - C_i
  // then we need the gradients \nabla r_i, and the approximate
  // Hessian is given by J'J, with J = dr_i/dx_j

  // we use a quasi Newton with backtracking and Wolfe condition check
  double c1 = 1e-4;
  double rho = 0.8;

  // the approximate Hessian
  DenseMatrix<double> M(nvar, nvar);

  // the quasi Newton step
  DenseVector<double> p(nvar);

  // the residual
  //DenseVector<double> residual;

  // the gradient of the objective function
  DenseVector<double> grad_f;

  residual_wz(a, c, u, keating, residual, grad);

  //residual = c_calc;
  //residual -= c_ref;
  double f_val = residual.l2_norm();
  f_val *= 0.5 * f_val;

  grad.get_transpose(M);
  M.vector_mult(grad_f, residual);
  M.right_multiply(grad);

  M.lu_solve(grad_f, p);

  int iter = 0;
  int max_it = 100;
  for ( ; (p.l1_norm() > 1e-3) && (iter < max_it); ++iter)
  {

    DenseVector<double> kt_old(keating);
    double t = 1;

    keating.add(-t, p);

    residual_wz(a, c, u, keating, residual, grad);

    //residual = c_calc;
    //residual -= c_ref;

    double new_f = residual.l2_norm();
    new_f *= 0.5 * new_f;
    for (int k = 0; (new_f > f_val + c1 * t * grad_f.dot(p)) && (k < 10); ++k)
    {
      t = rho * t;
      DenseVector<double> kt_old(keating);
      double t = 1;

      keating = kt_old;
      keating.add(-t, p);

      residual_wz(a, c, u, keating, residual, grad);
      //residual = c_calc;
      //residual -= c_ref;

      new_f = residual.l2_norm();
      new_f *= 0.5 * new_f;
    }

    grad.get_transpose(M);
    M.vector_mult(grad_f, residual);
    M.right_multiply(grad);

    M.lu_solve(grad_f, p);

  }

  alpha_0() = keating(0);
  alpha_1() = keating(2 % nvar);
  beta_0() = keating(1);
  beta_1() = keating(3 % nvar);

  //cerr << "after " << iter << " iterations:\n";
  //cerr << c_calc << endl;

  cerr << keating << endl;
}


void
AutomaticKeating::residual_wz(double a, double c, double u,
    const DenseVector<double>& keating,
    DenseVector<double>& residual,
    DenseMatrix<double>& gradients)
{

  const int nvar = keating.size();

  typedef F<double> double_t;

  vector<double_t> x(nvar);

  for (int i = 0; i < nvar; i++)
  {
    x[i] = keating(i);
    x[i].diff(i, nvar);
  }

  vector<double_t> y;
  keating_wz(a, c, u, x, y);

  unsigned int N = y.size();
  gradients.resize(N, nvar);
  residual.resize(N);


  vector<double> c_ref = {_c11, _c12, _c33, _c13, _c44, 0.5 * (_c11 - _c12)};

  for (int i = 0; i < N; i++)
  {
    double w = _weights[i];
    y[i] *= w;
    residual(i) = y[i].x() - w * c_ref[i];
    for (int j = 0; j < nvar; j++)
      gradients(i, j) = y[i].d(j);
  }
}


void
AutomaticKeating::calculate_zb_alpha(void)
{
  double alpha = (_c11 + 3.0 * _c12) * (get_a()  / 4.0);
  alpha_0() = alpha;
  alpha_1() = alpha;
}


void
AutomaticKeating::calculate_zb_beta(void)
{

  double beta = (_c11 - _c12) *  (get_a()  / 4.0);
  beta_0() = beta;
  beta_1() = beta;
}




