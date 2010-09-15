// $Id$

#ifndef _DDINTERFACEMODEL_H_
#define _DDINTERFACEMODEL_H_

#include "PhysicalModel.h"

#include <cassert>
#include <set>


class DriftDiffusionProperties;
class RecombinationModelInterface;
class Trap;

/*!
 * \brief The base class for the Drift-Diffusion boundary models
 *
 * Boundary conditions for the Drift-Diffusio equations are generally
 * written as Robin boundary conditions
 * \f[\beta\nabla \vec{\Phi}\cdot\vec{n} = g(u) - a(u)u\f
 * where \f$u\f$ is a potential, \f$\vec{\Phi}\f$ is the flux driven by
 * the gradient of \f$u\f$ and \f$\vec{n}\f$ is the outer normal
 * on the surface. The \f$\beta \in \{0, 1\}\f$ in the above equation is used
 * to allow for Dirichlet boundary conditions. For \f$\beta = 1\f$, the boundary
 * condition is of Robin type, for \f$\beta = 0\f$, it is of Dirichlet type.
 */
class DDInterfaceModel : public PhysicalModel
{

  public:

    //! The type of boundary condition
    enum BCType
    {
      DIRICHLET = 0,    //!< Dirichlet boundary condition
      ROBIN     = 1,    //!< Robin (type 3) boundary condition
      NEUMANN   = 1     //!< Neumann (type 2) boundary condition
    };

    //! Destructor
    virtual ~DDInterfaceModel(void) {};


    //! Create an interface model
    static DDInterfaceModel* create(const ModelOptions& options);


    //! Set the drift-diffusion properties object
    void set_dd_properties(DriftDiffusionProperties* ddprop);


    //! Compute the coefficients and their derivatives
    void compute();


    /*!
     * \brief Get the coefficients a
     *
     * Returns the coefficients \f$a(u)\f$ computed in compute()
     * for each equation.
     */
    const std::vector<double>& get_a(void) const;


    /*!
     * \brief Get the coefficients g
     *
     * Returns the coefficients \f$g(u)\f$ computed in compute()
     * for each equation.
     */
    const std::vector<double>& get_g(void) const;


    //! Get the coefficients b (to know the type of boundary condition)
    const std::vector<BCType>& get_b(void) const;


    /*!
     * \brief Get the coefficient derivatives with respect to the potentials
     *
     * Returns the partial derivatives (Jacobian)
     *  \f$\partial (g(u) - a(u)u)/\partial u\f$.
     */
    const std::vector<std::vector<double> >& get_jacobian(void) const;


    /*!
     * \brief Jacobian row for equation \c i
     */
    const std::vector<double>& get_jacobian_row(unsigned int i) const;


    //! Get the boundary condition type of variable \c var
    BCType get_type(unsigned int var) const;


    //! \c true if this is on an internal boundary
    bool is_internal_boundary(void) const;


    //! Tell it if it is an internal boundary
    void internal_bondary(bool is_internal_boundary);


    //! \c true if this boundary has non-zero current
    bool has_current(void) const;



  protected:

    //! Constructor
    DDInterfaceModel(const ModelOptions& options);


    /*!
     * \copydoc ElectricalContact::do_init()
     * Call this method from reimplementations!
     */
    virtual void do_init(void);


    //! Set the BC type for variable \c i
    void set_type(unsigned int var, BCType type);


    //! Do the actual calculation
    /*!
     * For Neumann/Robin boundary conditions, coeff_g and the jacobian
     * will be zeroed before calling do_compute().
     *
     */
    virtual void do_compute(void) {};


    //! Get the coefficients a
    double& coeff_a(unsigned int i);

    //! Get the coefficients b
    BCType& coeff_b(unsigned int i);

    //! Get the coefficients g
    double& coeff_g(unsigned int i);

    //! Get the jacobian
    double& jacobian(unsigned int i, unsigned int j);

    //! Tell the model that it has non-zero current
    void has_current(bool hascurrent);

    //! Get the DriftDiffusionProperties object
    DriftDiffusionProperties* get_dd_properties(void) const;


  private:


    //! The coefficients a
    std::vector<double> _coeff_a;

    //! The coefficients b
    std::vector<BCType> _coeff_b;

    //! The coefficients g
    std::vector<double> _coeff_g;

    //! The Jacobian (derivatives with respect to the potentials)
    std::vector<std::vector<double> > _jacobian;


    //! \c true if this is an internal boundary
    bool _internal_bd;


    //! \c true if this boundary has non-zero current
    bool _has_current;


    //! The driftdiffusion properties
    DriftDiffusionProperties* _ddprop;


    //! The electron traps
    std::set<Trap*> _etraps;

    //! The hole traps
    std::set<Trap*> _htraps;

    //! Recombination models
    std::set<RecombinationModelInterface*> _recombination_models;


    //! Calculate the trap contributions
    void _calculate_traps(double& q, double& dq_dEfn, double& dq_dEfp);


    //! calculate the recombinations
    void _calculate_recombination(double rec[6]);


    //! The creation method
    static TiberModelObject* _create(const ModelOptions& options)
    {
      return new DDInterfaceModel(options);
    }

    //! The destruction method
    static void _destroy(TiberModelObject* p)
    {
      delete p;
    }

};


//
// inline members
//

inline
void
DDInterfaceModel::set_dd_properties(DriftDiffusionProperties* ddprop)
{
  _ddprop = ddprop;
}


inline
const std::vector<double>&
DDInterfaceModel::get_a(void) const
{
  return _coeff_a;
}



inline
const std::vector<DDInterfaceModel::BCType>&
DDInterfaceModel::get_b(void) const
{
  return _coeff_b;
}



inline
const std::vector<double>&
DDInterfaceModel::get_g(void) const
{
  return _coeff_g;
}




inline
const std::vector<std::vector<double> >&
DDInterfaceModel::get_jacobian(void) const
{
  return _jacobian;
}


inline
const std::vector<double>&
DDInterfaceModel::get_jacobian_row(unsigned int i) const
{
  assert(i < _jacobian.size());
  return _jacobian[i];
}


inline
void
DDInterfaceModel::set_type(unsigned int var, BCType type)
{
  assert(var < _coeff_b.size());
  _coeff_b[var] = type;
  if (type == DIRICHLET)
  {
    _coeff_a[var] = 1.0;
    _jacobian[var][var] = -1;
  }
  else if (type == NEUMANN)
    _coeff_a[var] = 0.0;
}


inline
DDInterfaceModel::BCType
DDInterfaceModel::get_type(unsigned int var) const
{
  assert(var < _coeff_b.size());
  return _coeff_b[var];
}


inline
double&
DDInterfaceModel::coeff_a(unsigned int i)
{
  return _coeff_a[i];
}


inline
DDInterfaceModel::BCType&
DDInterfaceModel::coeff_b(unsigned int i)
{
  return _coeff_b[i];
}


inline
double&
DDInterfaceModel::coeff_g(unsigned int i)
{
  return _coeff_g[i];
}


inline
double&
DDInterfaceModel::jacobian(unsigned int i, unsigned int j)
{
  return _jacobian[i][j];
}


inline
bool
DDInterfaceModel::has_current(void) const
{
  return _has_current;
}


inline
void
DDInterfaceModel::has_current(bool hascurrent)
{
  _has_current = hascurrent;
}


inline
DriftDiffusionProperties*
DDInterfaceModel::get_dd_properties(void) const
{
  return _ddprop;
}


inline
bool
DDInterfaceModel::is_internal_boundary(void) const
{
  return _internal_bd;
}


inline
void
DDInterfaceModel::internal_bondary(bool is_internal_boundary)
{
  _internal_bd = is_internal_boundary;
}



#endif // _DDINTERFACEMODEL_H_
