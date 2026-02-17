/*  
 * This file is part of the tiberCAD module driftdiffusion.
 *
 * tiberCAD modules are licensed under the GNU General Public License v3.
 *
 * tiberCAD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published
 * by the Free Software Foundation, either version 3 of the License,
 * or (at your option) any later version.
 *
 * tiberCAD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with tiberCAD. If not, see <http://www.gnu.org/licenses/>.
 */

/*!
 * \file DDInterfaceModel.h
 * \brief tiberCAD driftdiffusion module header.
 *
 * \note This file is part of module driftdiffusion.
 */


#ifndef TC_DDINTERFACEMODEL_H
#define TC_DDINTERFACEMODEL_H

#include "DriftDiffusionProperties.h"
#include "tibercad/module/SolutionProvider.h"

#include "point.h"

#include <cassert>
#include <set>
#include "tibercad/base/libMeshDefs.h"

class DDBulkModel;
class MaterialBoundary;
class SimulationInterface;
class RecombinationModelInterface;
class Trap;
class FowlerNordheim;

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
class DDInterfaceModel : public DriftDiffusionProperties
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
    static DDInterfaceModel* create(const MaterialBoundary* boundary,
        const ModelOptions& options);


    //! reinitialize for the current element side
    void reinit(const Elem* elem, int side);


    //! Set the current side number and face normal
    void set_face_normal(const Point& n);


    //! Set a reference fermi level
    void set_reference_fermi_potentials(double fermi_e, double fermi_h);


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


    //! \c true if field emission model is defined
    bool has_field_emission(void) const;


    //! Get the field emission current
    FowlerNordheim* get_field_emission_model(void);


    //! Get the electron flux simulation
    SimulationInterface* get_eflux_simulation(void) const;

    //! Get the electron flux simulation
    SimulationInterface* get_hflux_simulation(void) const;


  protected:

    //! Constructor
    DDInterfaceModel(const ModelOptions& options);


    /*!
     * \copydoc ElectricalContact::do_init()
     * Call this method from reimplementations!
     */
    virtual void do_init(void);


    //! Create some of the submodels
    //virtual void prepare_submodels(void);


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

    //! Get the bulk DriftDiffusionProperties object
    DriftDiffusionProperties* get_bulk_dd_properties(void) const;

    //! Get the current side number
    int get_side_number(void) const;

    //! Get the current face normal
    const Point& get_face_normal(void) const;

    void get_reference_fermi_potentials(double& fermi_e, double& fermi_h) const;


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


    //! The current side number
    int _side;

    //! The current face normal
    Point _normal;

    //! A reference electron Fermi potential
    double _ref_fermi_e;

    //! A reference hole Fermi potential
    double _ref_fermi_h;

    //! The electron traps
    std::set<Trap*> _etraps;

    //! The hole traps
    std::set<Trap*> _htraps;

    //! Recombination models
    //std::set<RecombinationModelInterface*> _recombination_models;

    //! If we use Fowler-Nordheim emission
    FowlerNordheim* _emission;


    //! The electron flux in \f$cm^{-2}\f$
    double _eflux;
   
    //! The hole flux in \f$cm^{-2}\f$
    double _hflux;

    //! Solution provider for eflux
    SolutionProvider _eflux_sim;

    //! Solution provider for hflux
    SolutionProvider _hflux_sim;
    
    //! True if flux controlled for electrons
    bool _eflux_controlled;

    //! True if flux controlled for electrons
    bool _hflux_controlled;

    //! True if a flux predictor should be used
    bool _flux_predictor;

    //! The DD properties for material A
    DDBulkModel* _ddprop_A;

    //! The DD properties for material B
    DDBulkModel* _ddprop_B;


    //! The creation method
    static TiberModelObject* _create(const ModelOptions& options, const void*)
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
DDInterfaceModel::set_face_normal(const Point& n)
{
  _normal = n;
}


inline
void
DDInterfaceModel::set_reference_fermi_potentials(double fermi_e, double fermi_h)
{
  _ref_fermi_e = fermi_e;
  _ref_fermi_h = fermi_h;
}

inline
void
DDInterfaceModel::get_reference_fermi_potentials(double& fermi_e,
    double& fermi_h) const
{
  fermi_e = _ref_fermi_e;
  fermi_h = _ref_fermi_h;
}



inline
int
DDInterfaceModel::get_side_number(void) const
{
  return _side;
}


inline
const Point&
DDInterfaceModel::get_face_normal(void) const
{
  return _normal;
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
bool
DDInterfaceModel::has_field_emission(void) const
{
  return (_emission != NULL);
}


inline
FowlerNordheim*
DDInterfaceModel::get_field_emission_model(void)
{
  return _emission;
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


inline
SimulationInterface*
DDInterfaceModel::get_eflux_simulation(void) const
{
  return(_eflux_sim.simulation());
}

inline
SimulationInterface*
DDInterfaceModel::get_hflux_simulation(void) const
{
  return(_hflux_sim.simulation());
}


#endif // TC_DDINTERFACEMODEL_H
