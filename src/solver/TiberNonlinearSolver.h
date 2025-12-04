/*
 * This file is part of tiberCAD.
 *
 * tiberCAD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * tiberCAD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with tiberCAD. If not, see <https://www.gnu.org/licenses/>.
 */

/*!
 * \file TiberNonlinearSolver.h
 * \brief Private tiberCAD header.
 *
 * \internal
 */



#ifndef _TIBERNONLINEARSOLVER_H_
#define _TIBERNONLINEARSOLVER_H_

#include "tibercad/base/libMeshDefs.h"

// Libmesh includes
#include "nonlinear_solver.h"

class ModelOptions;
class XMonitor;

USELIBMESHTYPE(NumericVector);
USELIBMESHTYPE(SparseMatrix);


//! The TiberCAD nonlinear solver interface
/*!
 * This class provides the TiberCAD interface 
 * nonlinear solvers.
 * It is derived from the libmesh NonlinearSolver class
 *
 */
class TiberNonlinearSolver : public libMesh::NonlinearSolver<double>
{

  public:

    typedef libMesh::NonlinearImplicitSystem sys_type;

    //!  Constructor. Initializes  data structures
    TiberNonlinearSolver(sys_type& s);


    //! Destructor.
    virtual ~TiberNonlinearSolver(void);

    
    //! Create a linear solver
    //static TiberNonlinearSolver* create(const ModelOptions& options);

    
    //! Set options
    /*!
     * Call this method before calling solve()
     * 
     * Unspecified options are set to there default values,
     * \em not to their current values!
     */
    void set_options(const ModelOptions& options);

    
    //! Solve the system
    virtual std::pair<unsigned int, Real> solve(
        SparseMatrix<double>& jacobian,
        libMesh::NumericVector<double>& solution,
        libMesh::NumericVector<double>& residual,
        const double rtol,
        const unsigned int iter)
    {
      return this->solve(jacobian, solution, residual);
      static_cast<int>(rtol);
      static_cast<int>(iter);
    }

    
    //! Solve the system
    /*!
     * Use this method in TiberCAD instead of the libMesh defined ones
     */
    virtual std::pair<unsigned int, Real> solve(
        SparseMatrix<double>& jacobian,
        libMesh::NumericVector<double>& solution,
        libMesh::NumericVector<double>& residual) = 0;
    

    //! Get the relative nonlinear tolerance
    double get_nonlinear_rtol(void) const;

    //! Get the absolute nonlinear tolerance
    double get_nonlinear_atol(void) const;

    //! Get the absolute nonlinear tolerance
    double get_nonlinear_stol(void) const;

    //! Get the maximum number of nonlinear iterations
    int get_nonlinear_max_it(void) const;


    //! Set the relative nonlinear tolerance
    void set_nonlinear_rtol(double rtol);

    //! Set the absolute nonlinear tolerance
    void set_nonlinear_atol(double atol);

    //! Set the absolute nonlinear tolerance
    void set_nonlinear_stol(double stol);

    //! Set the maximum number of nonlinear iterations
    void set_nonlinear_max_it(int max_it);


    //! Set a pointer to the X monitor to be used
    void set_xmonitor(XMonitor* xmonitor);

    //! Get the X monitor
    XMonitor* get_xmonitor(void);


    //! Add a point to the X monitor
    /*!
     * \param iteration the iteration number
     * \param err the error
     * \param logarithm if \c true, plot \c log10(error)
     */
    void draw_point(double iteration, double error, bool logarithm = true);



    
  protected:

    //! Parse the options for solver specific stuff
    virtual void parse_options(const ModelOptions& options);

    


  private:

    //! The relative linear tolerance
    double _nonlinear_rtol;

    //! The absolute linear tolerance
    double _nonlinear_atol;

    //! The step tolerance
    double _nonlinear_stol;
      
    //! The maximum number of nonlinear iterations
    int _nonlinear_max_it;

    //! The X monitor
    XMonitor* _xmonitor;


};



//
// inline members
//



inline
TiberNonlinearSolver::~TiberNonlinearSolver(void)
{
}


inline
void
TiberNonlinearSolver::parse_options(const ModelOptions&)
{
}


inline
double
TiberNonlinearSolver::get_nonlinear_rtol(void) const
{
  return _nonlinear_rtol;
}


inline
double
TiberNonlinearSolver::get_nonlinear_atol(void) const
{
  return _nonlinear_atol;
}


inline
double
TiberNonlinearSolver::get_nonlinear_stol(void) const
{
  return _nonlinear_stol;
}



inline
int
TiberNonlinearSolver::get_nonlinear_max_it(void) const
{
  return _nonlinear_max_it;
}




inline
void 
TiberNonlinearSolver::set_nonlinear_rtol(double rtol)
{
  _nonlinear_rtol = rtol;
}

inline
void
TiberNonlinearSolver::set_nonlinear_atol(double atol)
{
  _nonlinear_atol = atol;
}

inline
void
TiberNonlinearSolver::set_nonlinear_stol(double stol)
{
  _nonlinear_stol = stol;
}

inline
void
TiberNonlinearSolver::set_nonlinear_max_it(int max_it)
{
  _nonlinear_max_it = max_it;
}




inline
XMonitor*
TiberNonlinearSolver::get_xmonitor(void)
{
  return _xmonitor;
}


inline
void
TiberNonlinearSolver::set_xmonitor(XMonitor* xmonitor)
{
  _xmonitor = xmonitor;
}


#endif // _TIBERNONLINEARSOLVER_H_
