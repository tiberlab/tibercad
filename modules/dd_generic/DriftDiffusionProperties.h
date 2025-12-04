/*  
 * This file is part of the tiberCAD module dd_generic.
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
 * \file DriftDiffusionProperties.h
 * \brief tiberCAD dd_generic module header.
 *
 * \note This file is part of module dd_generic.
 */


#ifndef _DRIFTDIFFUSIONPROPERTIES_H_
#define _DRIFTDIFFUSIONPROPERTIES_H_



#include "tibercad/physics/PhysicalModel.h"

#include "tibercad/physics/TemperatureInterface.h"
#include "CarrierProperties.h"
#include "tibercad/base/SimulationOptions.h"
#include "DriftDiffusionDefs.h"
#include "tibercad/base/TiberCad.h"
#include "tibercad/physics/Constants.h"
#include "tibercad/base/TypeDefs.h"
#include "tibercad/math/Tensor2.h"
#include "tibercad/base/RuntimeException.h"


#include "libmesh/tensor_value.h"
#include "libmesh/vector_value.h"
#include "libmesh/point.h"


#include <vector>
#include <set>
#include <map>


// forward declarations
class Dopant;
class Trap;
class SimulationInterface;
class DriftDiffusion;
class RecombinationModelInterface;
class ParticleDensity;


//! The base class for all drift-diffusion related semiconductor models
class DriftDiffusionProperties : public PhysicalModel
{

  private:

    //! An iterator for the recombination models
    typedef std::multimap<ID, RecombinationModelInterface*>::iterator
      recomb_iterator;

    //! A const iterator for the recombination models
    typedef std::multimap<ID, RecombinationModelInterface*>::const_iterator
      const_recomb_iterator;

  public:

    class RecombinationModelIterator : public recomb_iterator
    {
      public:
        RecombinationModelIterator(const recomb_iterator& it) :
          recomb_iterator(it) {}

        RecombinationModelIterator(const RecombinationModelIterator& it) :
          recomb_iterator(it) {}

        ID id(void) { return((*this)->first); }

        RecombinationModelInterface* operator*(void) { return((*this)->second); }

    };

    //! A nested class that holds all point data
    class PointData
    {
      public:

        //! Constructor
        PointData(void);

        //! The electric potential
        double electric_potential = 0;


        //! The electric potential of the step before
        double old_electric_potential = 0;


        //! The ionized donor density
        double ionized_donor_density = 0;

        //! The ionized donor density derivative
        double ionized_donor_density_derivative = 0;

        //! The ionized acceptor density
        double ionized_acceptor_density = 0;

        //! The ionized acceptor density derivative
        double ionized_acceptor_density_derivative = 0;


        //! The trapped electron density
        double ionized_electron_traps = 0;

        //! The trapped hole density
        double ionized_hole_traps = 0;

        //! The trapped density derivatives d/dfermi_e, d/dfermi_h
        std::vector<double> ionized_traps_derivative;

        //! The total charge density
        double charge_density = 0;


        //! A map storing the fermi potential for each carrier identified by id
        std::vector<double> fermi_potential;

        //! Old fermi potential
        std::vector<double> old_fermi_potential;

        //! A map storing carrier densities
        std::vector<double> q_density;

         //! A map storing carrier density derivatives
        std::vector<double> q_density_derivative;

        //! Total charge density derivatives
        std::vector<double> charge_density_derivative;

        //! A map storing carrier conductivities
        std::vector<double> q_conductivity;

        //! A map storing carrier conductivities
        std::vector<double> q_conductivity_derivative_qFermi;

        //! A map storing carrier mobilities
        std::vector<double> q_mobility;

        //! A map storing net recombination rates
        std::vector<double> q_recombination_rate;

        //! The net recombination rate derivatives w.r.t potentials
        /*!
         * order is the same as the variable order
         */
        std::vector<std::vector<double>> q_recombination_rate_derivatives;

        //! The carrier temperature in eV (\f$= k_B T_e / e\f$)
        std::vector<double> carrier_vt;

    };


    //! An ID to identify unknown carriers
    static const ID unknown_carrier_id = static_cast<ID>(-1);



    //! A default (empty) destructor.
    virtual ~DriftDiffusionProperties(void);


    //! Create a named drift-diffusion model
    /*!
     * The model is created according to the given model name.
     * If it is not known, the NULL pointer is returned.
     *
     * \param name the model name
     * \param options the options as given in the input file
     * \return a pointer to the newly created object
     */
    static DriftDiffusionProperties* create(const std::string& name,
        const Material* mat,
        const ModelOptions& options = ModelOptions());


    //! Set the list of known carriers, defining their ordering
    /*!
     * Carrier properties could be read from the input in any order, however the module
     * will assume a particular order to address them. This is assured by passing all
     * known carrier names in the same order as they are used in \c DriftDiffusion,
     * immediately after creating the model.
     *
     * \param known_carriers the sequence of carrier names in global ordering
     */
    void set_known_carriers(const std::vector<std::string>& known_carriers);


    //! Get the ordered vector of known carriers
    const std::vector<std::string>& get_known_carriers(void) const;

    //! Get the number of known carriers
    unsigned int n_known_carriers(void) const;


    //! Get the ID of a carrier
    /*!
     * \return \c unknown_carrier_id if \c name is not present
     *
     * This method should be called only with carrier names known
     * to be present.
     */
    int get_carrier_id(const std::string& name) const;


    // ! Lock the parameters
    /* !
     * Call this method before reinit()
     *
     * You can provide a pointer to a PointData object created outside
     * or let it create one internally
     */
    //void lock(PointData* pd = NULL);


    // ! Unlock the parameters
    /* !
     * Call this after all calculations on the element hav been done
     */
    //void unlock(void);


    // ! (Re-)Initialize for the given element
    /* !
     * \c reinit() calls \c prepare_element_data() which needs to be
     * implemented in derived classes
     */
    //void reinit(const Elem* elem);


    //! Set the polarization vector
    //void set_polarization(const RealVectorValue& polarization);


    //! Set the coordinates
    void set_coordinates(const Point& p);

    //! Set the carrier temperatures
    /*!
     * The carrier temperatures have to be set before
     * calling any of the \c calculate_xxx methods
     * This method sets thesame temperature for all carriers.
     *
     * \param T the carrier temperature
     */
    void set_carrier_temperatures(double T);

    //! Set the temperature of carrier \c id
    void set_carrier_temperature(ID carrier, double T);


    //! Set the electric field
    /*!
     * The electric field has to be set before calling any of the
     * \c calculate_xxx methods.
     *
     * \param E the electric field
     */
    void set_electric_field(const libMesh::RealGradient& E);
    void set_old_electric_field(const libMesh::RealGradient& E);


    //! Get the electric field
    const libMesh::RealGradient& get_electric_field(void) const;
    const libMesh::RealGradient& get_old_electric_field(void) const;


    //! Get the element we are currently working on
    const Elem* get_element(void) const;

    //! Get the coordinates of the point we are currently working on
    const Point& get_coordinates(void) const;

    //! Set the strain
    void set_strain(const Tensor2& strain);

    //! Get the strain
    const Tensor2& get_strain(void) const;


    //! Set the lattice temperature (in K)
    void set_lattice_temperature(double T);

    //! Set the lattice temperature (in eV)
    void set_lattice_vt(double vt);

    //! Get the lattice temperature (in units of eV)
    double get_lattice_temperature(void) const;



    //! Calculate electron and hole densities and derivatives
    /*!
     * This method calculates electron and hole densities and their
     * derivatives with respect to the electric and the
     * electro-chemical potentials.
     */
    void calculate_densities(void);


    //! Calculated trapped particle densities
    void calculate_traps(void);


    //! Calculate the ionized dopant densities and derivatives
    /*!
     * The total density of ionized donors and acceptors is calculated,
     * respectively. Their derivative with respect to the electric
     * potential are also computed (the derivative with respect to the
     * electro-chemical potentials have the same value but opposite sign).
     */
    void calculate_ionized_dopants(void);


    //! Calculate net recombination rates and derivatives
    /*!
     * The recombination models calculate the derivatives with respect to
     * the densities.
     */
    void calculate_net_recombination_rates(void);



    //! Get the ionized donor density
    /*!
     * \return the ionized donor density
     */
    double get_ionized_donor_density(void) const
      { return _pd->ionized_donor_density; };


    //! Get the ionized donor density derivative
    /*!
     * \return the ionized donor density derivative with respect to the
     * electric potential
     */
    double get_ionized_donor_density_derivative(void) const
      { return _pd->ionized_donor_density_derivative; };


    //! Get the ionized acceptor density
    /*!
     * \return the ionized acceptor density
     */
    double get_ionized_acceptor_density(void) const
      { return _pd->ionized_acceptor_density; };


    //! Get the ionized acceptor density derivative
    /*!
     * \return the ionized acceptor density derivative with respect to the
     * electric potential
     */
    double get_ionized_acceptor_density_derivative(void) const
      { return _pd->ionized_acceptor_density_derivative; };


    //! Get the ionized electron trap density
    /*!
     * \return the ionized electron trap density
     */
    double get_ionized_electron_traps(void) const
      { return _pd->ionized_electron_traps; };


    //! Get the ionized hole trap density
    /*!
     * \return the ionized hole trap density
     */
    double get_ionized_hole_traps(void) const
      { return _pd->ionized_hole_traps; };


    //! Get the total charge density
    /*!
     *
     * \f$ \rho = p - n + N_D^+ - N_A^- \f$
     *
     * \return the total charge density \f$\rho\f$
     *
     * \note
     * The charge density is returned in units of the elementary charge
     * \f$e\f$, \em not in Coulomb!
     */
    double get_charge_density(void) const;


    //! Get the derivatives of the total charge density
    /*!
     * The derivatives are given w.r.t. the quasi Fermi levels
     */
    void get_charge_density_derivatives(std::vector<double>& derivatives) const;

    //! Get the derivative of the total charge density
    /*!
     * The derivative is given w.r.t. the quasi Fermi level
     *
     * \param id the id of the particle
     */
    double get_charge_density_derivative(ID id) const;


    //! Get the total electric polarization
    /*!
     * The total electric polarization \b P is the sum of the
     * pyroelectric and piezoelectric polarization
     */
    //const RealVectorValue& get_total_polarization(void) const
    //  { return _polarization; };


    //! Get the relative permittivity tensor
    const libMesh::RealTensor& get_relative_permittivity(void) const
      { return _permittivity; };


    //! Get equilibrium fermi level
    double get_equilibrium_fermi_level(void) const
      { return _equilibrium_fermi_level; };


    //void get_net_recombination_rates(std::vector<double>& rates);
    //! Get the iterator to the first recombination model
    RecombinationModelIterator recombination_models_begin(void);


    //! Get the past-the-end iterator for the recombination models
    RecombinationModelIterator recombination_models_end(void);

    //! Get the IDs of the registered recombination models
    /*!
     * \return number of registered models
     */
    int get_net_recombination_rate_IDs(std::vector<ID>& ids);


    //! Get the recombination model with ID \c id
    RecombinationModelInterface* get_recombination_model(ID id);


    //! Get the recombination model with ID \c id
    const RecombinationModelInterface* get_recombination_model(ID id) const;


    //! get the net recombination rate of model \c id
    std::pair<double, double> get_net_recombination_rate(ID id);


    //! Clear all recombination rates
    void clear_recombination(void);


    //! Returns the number of recombination models
    int get_number_of_recombination_models(void) const;


    //! Get the electric potential
    double get_electric_potential(void) const;

    double get_old_phi(void) const
      { return _pd->old_electric_potential; };


    //! Get the point data
    const PointData& get_point_data(void) const
        { return get_pd(); }


    //! Create the recombination models
    void create_recombination_models(void);


    //! Tells if we are doing equilibrium calculation
    bool has_solution(void) const;


    //! Set the potentials
    void set_potentials(double phi, const std::vector<double> & q_phif = std::vector<double>(0));
    void set_el_potential(double phi);
    void set_fermi_potential(ID id, double phif);

    //! Set the old potentials
    void set_old_potentials(double phi, const std::vector<double> & q_phif);

    void set_old_el_potential(double phi);
    void set_old_fermi_potential(ID id, double phif);

    //! Set the gradient electro-chemical potentials
    void set_grad_fermi(const std::vector<libMesh::RealGradient> & q_gradf);
    void set_grad_fermi(ID id, const libMesh::RealGradient & gradf);

    //! Get the gradient of the electrochemical potential
    const std::vector<libMesh::RealGradient>& get_grad_fermi(void) const;
    void get_grad_fermi(ID id, libMesh::RealGradient& gradf) const;

    //! Get carrier densities
    const std::vector<double>& get_q_density(void) const
      { return _pd->q_density; };
    double get_q_density(ID id) const;

    //! Get carrier density derivatives
    const std::vector<double>& get_q_density_derivative(void) const
      { return _pd->q_density_derivative; };
    double get_q_density_derivative(ID id) const;

    //! Get carrier conductivities
    const std::vector<double>& get_q_conductivity(void) const
      { return _pd->q_conductivity; };
    double get_q_conductivity(ID id) const;

    //! Get carrier mobilities
    const std::vector<double>& get_q_mobility(void) const
      { return _pd->q_mobility; };
    double get_q_mobility(ID id) const;


    //! Get the net recombination rates
    /*!
     * Get \f$R_{net} = R - G\f$ as
     * calculated by \c calculate_all(...)
     *
     */
    const std::vector<double>&  get_net_q_recombination_rate(void) const
      { return _pd->q_recombination_rate; };
    double get_net_q_recombination_rate(ID id) const;

    //! Get the net recombination rate derivatives
    const std::vector<std::vector<double>>&
            get_net_q_recombination_rate_derivatives(void) const
              { return _pd->q_recombination_rate_derivatives; };
    const std::vector<double>& get_net_q_recombination_rate_derivatives(ID id) const;


    //! Get the carrier electro-chemical potential
    const std::vector<double>& get_q_fermi_potential(void) const
      { return _pd->fermi_potential; };
    double get_q_fermi_potential(ID id) const;

    //! Get the carrier electro-chemical potential
    const std::vector<double>& get_old_q_fermi_potential(void) const
      { return _pd->old_fermi_potential; };
    double get_old_q_fermi_potential(ID id) const;

    //! Get carrier properties
    /*!
     * This method will only return carrier properties present in the object.
     */
    const std::map<ID, CarrierProperties*>& get_carrier_properties(void) const
      { return _carrier_properties; };

    //! Get carrier properties
    /*!
     * This method will only return carrier properties present in the object.
     */
    std::map<ID, CarrierProperties*>& get_carrier_properties(void)
      { return _carrier_properties; };

    const CarrierProperties* get_carrier_properties(ID id) const;

    CarrierProperties* get_carrier_properties(ID id);

    CarrierProperties* get_carrier_properties(const std::string& name) const;

    //! Get the carrier band edge
    double get_carrier_band_edge(ID id) const;



  protected:



    //! The empty constructor.
    DriftDiffusionProperties(const ModelOptions& options);


    //! Initialize this model
    /*!
     * This reads the database and calls init for all submodels
     * A derived class which reimplements this method has to call
     * explicitly the one of this class!
     */
    virtual void do_init(void) override;

    //! Copy some stuff from the original model
    virtual void copy_from(const PhysicalModel* rhs) override;


    //! Create some of the submodels
    virtual void prepare_submodels(void) override;


    //! Set the element we are currently working on
    void set_element(const Elem* elem);


    //! This method gets called from reinit()
    /*!
     * It can be used to setup data that is constant in an element, e.g.
     * strain related stuff, band edges.
     * This method can be overriden by derived classes.
     */
    virtual void prepare_element_data(void) {};



    //! Get the point data structure
    PointData& get_pd(void) const;


    //! Get the strain as writable reference
    Tensor2& get_strain(void);


    //! Get the temperature interface
    TemperatureInterface& get_temperature_interface(void);

    //! Get the strain interface
    //StrainInterface& get_strain_interface(void);


    //! Tells if we should assume inhomogeneous band parameters
    bool is_inhomogeneous(void) const;






    //! Set the relative permittivity tensor
    void set_relative_permittivity(const libMesh::RealTensor& eps)
      { _permittivity = eps; };





    //! Set the equilibrium Fermi level
    void set_equilibrium_fermi_level(double Ef)
      { _equilibrium_fermi_level = Ef; }


    //! Set the bulk equilibrium density for carrier \c id
    void set_bulk_equilibrium_density(ID id, double density);

    //! Set carrier properties
    /*!
     * This is used e.g. from interface models.
     */
    void set_carrier_properties(const std::map<ID, CarrierProperties*>& cp);



  private:

    //! The interface to the lattice temperature simulation
    TemperatureInterface _lattice_temp;



    //! \c true if we should assume inhomogeneous band parameters
    bool _is_inhomogeneous;





    //! The equilibrium fermi level
    /*!
     * The fermi level such that \f$n=n_0,\,p=p_0\f$
     */
    double _equilibrium_fermi_level;


    //! The point-wise data
    PointData* _pd;



    //! The electric field
    libMesh::RealGradient _electric_field;
    libMesh::RealGradient _old_electric_field;


    // ! The total electric polarization
    //RealVectorValue _polarization;

    //! The relative permittivity tensor
    libMesh::RealTensor _permittivity;


    //! The copy constructor is disabled
    DriftDiffusionProperties(const DriftDiffusionProperties& rhs);


    //! The assignment operator is disabled
    DriftDiffusionProperties& operator=(const DriftDiffusionProperties& rhs);


    //! Add a recombination model
    /*!
     * Creates and adds a new recombination model from the given name and
     * options.
     *
     * \param model_name the name of the model
     * \param options the options for the model
     */
     void add_recombination_model(const std::string& model_name,
        const ModelOptions& options = ModelOptions());



    //! The element we are currently working on
    const Elem* _elem;

    //! The coordinates of the point we are working on
    Point _coord;


    //! The strain
    Tensor2 _strain;


    //! The electron traps
    std::set<Trap*> _etraps;

    //! The hole traps
    std::set<Trap*> _htraps;


    //! The recombination models
    std::multimap<ID, RecombinationModelInterface*> _recombination_models;


    // ! The polarization models
    //std::vector<PolarizationModel*> _pm;


    //! The lattice temperature in eV (\f$= k_B T_{lat} / e\f$)
    double _lattice_vt;


    //! The gradient of the electrochemical-potential
    std::vector<libMesh::RealGradient> _grad_fermi;

    //! The densities in equilibrium
    std::vector<double> _equilibrium_densities;

    //! The carrier properties
    /*!
     * Carrier properties are assumed to be elemental data, \em not nodal data
     */
    std::map<ID, CarrierProperties*> _carrier_properties;


    //! This vector contains all the known carriers in a consistent ordering
    /*!
     * The order is set from \c DriftDiffusion, and used also to define the
     * system and solution variables
     */
    std::vector<std::string> _known_carriers_ordered;

    //! The carrier relevant for donor properties
    unsigned int _donor_reference_carrier;

    //! The carrier relevant for acceptor properties
    unsigned int _acceptor_reference_carrier;

};


//
// inline members
//




inline
void
DriftDiffusionProperties::set_known_carriers(
    const std::vector<std::string>& known_carriers)
{
  _known_carriers_ordered = known_carriers;
}


inline
const std::vector<std::string>&
DriftDiffusionProperties::get_known_carriers(void) const
{
  return(_known_carriers_ordered);
}

inline
unsigned int
DriftDiffusionProperties::n_known_carriers(void) const
{
  return(_known_carriers_ordered.size());
}


inline
void
DriftDiffusionProperties::set_coordinates(const Point& p)
{
  _coord = p;
}



inline
void
DriftDiffusionProperties::set_bulk_equilibrium_density(ID id, double density)
{
  _equilibrium_densities[id] = density;
}


inline
double
DriftDiffusionProperties::get_electric_potential(void) const
{
  return _pd->electric_potential;
}



inline
void
DriftDiffusionProperties::set_electric_field(const libMesh::RealGradient& E)
{
  _electric_field = E;
}




inline
void
DriftDiffusionProperties::set_old_electric_field(const libMesh::RealGradient& E)
{
  _old_electric_field = E;
}



inline
const libMesh::RealGradient&
DriftDiffusionProperties::get_electric_field(void) const
{
  return _electric_field;
}

inline
const libMesh::RealGradient&
DriftDiffusionProperties::get_old_electric_field(void) const
{
  return _old_electric_field;
}



/*
inline
void
DriftDiffusionProperties::set_densities(double n, double p)
{
  _pd->electron_density = n;
  _pd->hole_density = p;
  _pd->electron_density_derivative = n / _pd->electron_vt;
  _pd->hole_density = p;
  _pd->hole_density_derivative = -p / _pd->hole_vt;
}
*/


inline
void
DriftDiffusionProperties::set_lattice_temperature(double T)
{
  _lattice_vt = T * Constants::k_B;
}


inline
void
DriftDiffusionProperties::set_lattice_vt(double vt)
{
  _lattice_vt = vt;
}


inline
double
DriftDiffusionProperties::get_lattice_temperature(void) const
{
  return _lattice_vt;
}

inline
void
DriftDiffusionProperties::set_element(const Elem* elem)
{
  _elem = elem;
}

inline
const Elem*
DriftDiffusionProperties::get_element(void) const
{
  return _elem;
}

inline
const Point&
DriftDiffusionProperties::get_coordinates(void) const
{
  return _coord;
}


inline
void
DriftDiffusionProperties::get_charge_density_derivatives(
    std::vector<double>& derivatives) const
{
  derivatives.resize(0);
  derivatives.resize(_carrier_properties.size(), 0.0);

  for (auto& cp: _carrier_properties)
    derivatives[cp.first] += _pd->charge_density_derivative[cp.first];

  //derivatives[0] = get_electron_density_derivative()
  //  - get_ionized_donor_density_derivative()
  //  - _pd->ionized_electron_traps_derivative;
  //derivatives[1] = -get_hole_density_derivative()
  //  + get_ionized_acceptor_density_derivative()
  //  - _pd->ionized_hole_traps_derivative;
}

inline
RecombinationModelInterface*
DriftDiffusionProperties::get_recombination_model(ID id)
{
  RecombinationModelInterface* rec = NULL;
  recomb_iterator it = _recombination_models.find(id);
  if (it != _recombination_models.end())
    rec = it->second;

  return rec;
}


inline
const RecombinationModelInterface*
DriftDiffusionProperties::get_recombination_model(ID id) const
{
  RecombinationModelInterface* rec = NULL;
  const_recomb_iterator it = _recombination_models.find(id);
  if (it != _recombination_models.end())
    rec = it->second;

  return rec;
}


inline
void
DriftDiffusionProperties::set_strain(const Tensor2& strain)
{
  _strain = strain;
}


inline
const Tensor2&
DriftDiffusionProperties::get_strain(void) const
{
  return _strain;
}


inline
Tensor2&
DriftDiffusionProperties::get_strain(void)
{
  return _strain;
}








//inline
//void
//DriftDiffusionProperties::set_polarization(const RealVectorValue& polarization)
//{
//  _polarization = _relax_polariz * polarization;
//}


inline
int
DriftDiffusionProperties::get_number_of_recombination_models(void) const
{
  return _recombination_models.size();
}

/*
inline
double
DriftDiffusionProperties::get_electron_thermoelectric_power(void) const
{
  return _eTEpower;
}


inline
double
DriftDiffusionProperties::get_hole_thermoelectric_power(void) const
{
  return _hTEpower;
}

inline
libMesh::RealGradient
DriftDiffusionProperties::get_electron_thermoelectric_power_gradient(void) const
{

  return  _eTEpowerGrad;

}

inline
libMesh::RealGradient
DriftDiffusionProperties::get_hole_thermoelectric_power_gradient(void) const
{

  return  _hTEpowerGrad;

}
*/


inline
DriftDiffusionProperties::PointData&
DriftDiffusionProperties::get_pd(void) const
{
  return *_pd;
}


inline
TemperatureInterface&
DriftDiffusionProperties::get_temperature_interface(void)
{
  return _lattice_temp;
}

//inline
//StrainInterface&
//DriftDiffusionProperties::get_strain_interface(void)
//{
//  return _strain_if;
//}


inline
bool
DriftDiffusionProperties::is_inhomogeneous(void) const
{
  return _is_inhomogeneous;
}


inline
DriftDiffusionProperties::RecombinationModelIterator
DriftDiffusionProperties::recombination_models_begin(void)
{
  return(RecombinationModelIterator(_recombination_models.begin()));
}


inline
DriftDiffusionProperties::RecombinationModelIterator
DriftDiffusionProperties::recombination_models_end(void)
{
  return(RecombinationModelIterator(_recombination_models.end()));
}

//New generic module

inline
void
DriftDiffusionProperties::set_potentials(double phi, const std::vector<double>& q_phif)
{
  _pd->electric_potential = phi;
  _pd->fermi_potential = q_phif;
  // if vector was empty:
  _pd->fermi_potential.resize(n_known_carriers(), 0.0);
}

inline
void
DriftDiffusionProperties::set_el_potential(double phi)
{
  _pd->electric_potential = phi;
}

inline
void
DriftDiffusionProperties::set_fermi_potential(ID id, double phif)
{
  _pd->fermi_potential.resize(n_known_carriers());
  _pd->fermi_potential[id] = phif;
}

inline
void
DriftDiffusionProperties::set_old_potentials(double phi, const std::vector<double> & q_phif)
{
  _pd->old_electric_potential = phi;
  _pd->old_fermi_potential.resize(n_known_carriers());
  _pd->old_fermi_potential = q_phif;
}

inline
void
DriftDiffusionProperties::set_old_el_potential(double phi)
{
  _pd->old_electric_potential = phi;
}

inline
void
DriftDiffusionProperties::set_old_fermi_potential(ID id, double phif)
{
  _pd->old_fermi_potential.resize(n_known_carriers());
  _pd->old_fermi_potential[id] = phif;
}

inline
void 
DriftDiffusionProperties::set_grad_fermi(const std::vector<libMesh::RealGradient> & q_gradf)
{
  _grad_fermi.resize(n_known_carriers());
  _grad_fermi = q_gradf;
}

inline
void
DriftDiffusionProperties::set_grad_fermi(ID id, const libMesh::RealGradient& gradf)
{
  _grad_fermi.resize(n_known_carriers());
  _grad_fermi[id] = gradf;

}

inline
const std::vector<libMesh::RealGradient>&
DriftDiffusionProperties::get_grad_fermi(void) const
{
  return _grad_fermi;
}

inline
void 
DriftDiffusionProperties::get_grad_fermi(ID id, libMesh::RealGradient& gradf) const
{
  gradf = _grad_fermi.at(id);
}

inline
double 
DriftDiffusionProperties::get_q_density(ID id) const
{
  return _pd->q_density[id];
}

inline
double 
DriftDiffusionProperties::get_q_density_derivative(ID id) const
{
  return _pd->q_density_derivative[id];
}


inline
double
DriftDiffusionProperties::get_charge_density_derivative(ID id) const
{
  return _pd->charge_density_derivative[id];
}

inline
double 
DriftDiffusionProperties::get_q_conductivity(ID id) const
{
  return _pd->q_conductivity[id];
}

inline
double 
DriftDiffusionProperties::get_q_mobility(ID id) const
{
  return _pd->q_mobility[id];
}

inline
double 
DriftDiffusionProperties::get_net_q_recombination_rate(ID id) const
{
  return _pd->q_recombination_rate[id];
}

inline
const std::vector<double>&
DriftDiffusionProperties::get_net_q_recombination_rate_derivatives(ID id) const
{
  return _pd->q_recombination_rate_derivatives[id];
}


inline
double 
DriftDiffusionProperties::get_q_fermi_potential(ID id) const
{
  return _pd->fermi_potential[id];
}

inline
double 
DriftDiffusionProperties::get_old_q_fermi_potential(ID id) const
{
  return _pd->old_fermi_potential[id];
}


inline
int
DriftDiffusionProperties::get_carrier_id(const std::string& name) const
{
  int id = unknown_carrier_id;
  for (unsigned int i = 0; i < _known_carriers_ordered.size(); ++i)
  {
    if (name == _known_carriers_ordered[i])
      id = i;
  }
  return(id);
}



inline
CarrierProperties*
DriftDiffusionProperties::get_carrier_properties(const std::string& name) const
{
  CarrierProperties* prop = nullptr;
  for (auto&& it : _carrier_properties)
  {
    if (name == it.second->get_name())
    {
      prop = it.second;
    }
  }
  return(prop);
}


inline
const CarrierProperties* 
DriftDiffusionProperties::get_carrier_properties(ID id) const
{
  const CarrierProperties* ptr = nullptr;
  auto it = _carrier_properties.find(id);
  if (it != _carrier_properties.end())
    ptr = it->second;

  return(ptr);
}

inline
CarrierProperties* 
DriftDiffusionProperties::get_carrier_properties(ID id)
{
  CarrierProperties* ptr = nullptr;
  auto it = _carrier_properties.find(id);
  if (it != _carrier_properties.end())
    ptr = it->second;

  return(ptr);
}

inline
double 
DriftDiffusionProperties::get_carrier_band_edge(ID id) const
{
  double energy = 0;
  auto it = _carrier_properties.find(id);
  if (it != _carrier_properties.end())
    energy = (it->second)->get_band_edge();

  return(energy);
}


inline
void 
DriftDiffusionProperties::set_carrier_temperatures(double T)
{
  _pd->carrier_vt.resize(_known_carriers_ordered.size());
  for (auto& cp : _pd->carrier_vt)
     cp = T;
}


inline
void 
DriftDiffusionProperties::set_carrier_temperature(ID id, double T)
{
  _pd->carrier_vt.resize(_known_carriers_ordered.size());
  _pd->carrier_vt[id] = T;
}



#endif /* _DRIFTDIFFUSIONPROPERTIES_H_ */

















