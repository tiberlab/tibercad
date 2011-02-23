// $Id$

#ifndef _EXCITONPROPERTIES_H_
#define _EXCITONPROPERTIES_H_


#ifndef TIBER_MODULE_PREFIX
# define TIBER_MODULE_PREFIX ex
#endif


#include "PhysicalModel.h"
#include "SimulationOptions.h"
#include "TemperatureInterface.h"
#include "TiberCad.h"

#include "vector_value.h"

#include <vector>

// forward declarations
class Point;
class Elem;

class TBDLEXPORT ExcitonProperties : public PhysicalModel
{

  public:

    //! A default (empty) destructor.
    virtual ~ExcitonProperties(void);


    //! Create a named model
    /*!
     * The model is created according to the given model name.
     * If it is not known, the NULL pointer is returned.
     *
     * \param name the model name
     * \param options the options as given in the input file
     * \return a pointer to the newly created object
     */
    static ExcitonProperties* create(const std::string& name,
        const ModelOptions& options = ModelOptions());

     //! Returns the exciton energy \f$ Eg - R \f$
     double get_exciton_energy(void);

    //! Set the statistics to be used
    /*!
     * \param statistics the statistics
     */
    void set_statistics(TiberCad::Statistics statistics);

    //! Get the statistics to be used
    /*!
     * \return the statistics
     */
    TiberCad::Statistics get_statistics(void) const;

    //! Set the mobility model
    /*!
     * Creates a mobility model for the excitons from the given model
     * name and options.
     */
    void set_mobility_model(const std::string& model_name,
        const ModelOptions& options = ModelOptions()) {};


    //! (Re-)Initialize for the given element
    /*!
     * \c reinit() calls \c prepare_element_data() which needs to be
     * implemented in derived classes
     *
     * \param elem the current element
     * \param dd_prop a pointer to the semiconductor model
     */
    void reinit(const Elem* elem);


    //! Set the coordinates
    void set_coordinates(const Point& p);

    //! Set the carrier temperature
    void set_carrier_temperature(double T_x);

    //! Set the effective exciton potential
    void set_effective_potential(double eff_potential);

    //! Set the exciton energy \f$E_x\f$
    /*!
     * \f$E_x=E_g-R\f$ where \f$E_g\f$ and \f$R\f$ are the band gap and
     * the Exciton binding energy, respectively.
     */
    void set_energy(double energy);

    //! Set the exciton DOS
    void set_density_of_states(double DOS);

    //! Set the carrier density
    void set_density(double x);


    //! Get the element we are currently working on
    const Elem* get_element(void) const;

    //! Get the coordinates of the point we are currently working on
    const Point& get_coordinates(void) const;


    //! Calculate the exciton density
    void calculate_density(void);

    //! Calculate the exciton net recombination rate
    void calculate_net_recombination_rate(void);

    //! Calculate the exciton mobility
    void calculate_mobility(void);


    //! Get the lattice temperature
    double get_lattice_temperature(void);

    //! Get the carrier temperature
    double get_carrier_temperature(void);

    //! Get the exciton density
    /*!
     * Get the exciton density as calculated by \c calculate_all(...)
     *
     * \return the exciton density
     */
    double get_density(void) const
      { return density; };

    //! Get the exciton density derivative
    /*!
     * \return the exciton density derivative with respect to the
     * electro-chemical potential
     */
    double get_density_derivative(void) const
      { return density_derivative; };


    //! Get the exciton recombination rate
    double get_net_recombination_rate(void) const;

    //! Get the exciton recombination rate derivative
    /*!
     * Get \f$\frac{\partial R}{\partial\phi_x}\f$
     */
    double get_net_recombination_rate_derivative(void) const
        { return recombination_rate_derivative; };

    //! Get the nonradiative recombination rate
    virtual double get_nonradiative_recombination_rate(void);

    //! Get the radiative recombination rate
    virtual double get_radiative_recombination_rate(void);

    //! Get the dissociation rate
    virtual double get_dissociation_rate(void);

    //! Get the exciton mobility
    /*!
     * \return the exciton mobility
     */
    double get_mobility(void) const
      { return mobility; };

  protected:

    //! The empty constructor.
    ExcitonProperties(const ModelOptions& options);


    //! Initialize this model
    /*!
     * This reads the database and calls init for all submodels
     * A derived class which reimplements this method has to call
     * explicitly the one of this class!
     */
    virtual void do_init(void) {};

    virtual void do_recombination(void) {};

    virtual void do_mobility(void) {};


    //! This method gets called from reinit()
    /*!
     * It can be used to setup data that is constant in an element, e.g.
     * strain related stuff, band edges.
     */
    virtual void prepare_element_data(void) {};


    //! The lattice temperature in eV (\f$= k_B T_{lat} / e\f$)
    double lattice_vt;

    //! The exciton temperature in eV (\f$= k_B T_h / e\f$)
    double exciton_vt;

    //! The effective exciton potential
    double eff_pot;

    //! The density
    double density;

    //! The density derivative
    double density_derivative;

    //! The mobility
    double mobility;

    //! The radiative recombination rate
    double rad_recomb_rate;

    //! The nonradiative recombination rate
    double nonrad_recomb_rate;

    //! The dissociation rate
    double dissoc_rate;

    //! The net recombination rate
    double net_recomb_rate;

    //! The derivative of the net recombination rate
    double recombination_rate_derivative;


  private:

    //! The copy constructor is disabled
    ExcitonProperties(const ExcitonProperties& rhs);

    //! The assignment operator is disabled
    ExcitonProperties& operator=(const ExcitonProperties& rhs);

    //! The interface to the lattice temperature simulation
    TemperatureInterface _lattice_temp;

    //! The element we are currently working on
    const Elem* _elem;

    //! The coordinates of the point we are working on
    const Point* _coord;

    //! The exciton energy
    double _energy;

    //! The DOS of the excitons
    double _DOS;

    //! The statistics to be used
    TiberCad::Statistics _statistics;

    //! Calculate the density and its derivative
    void calculate_density_and_derivative(double arg, double& density,
        double& derivative) const;

    //! Calculate the density for a given argument
    double calculate_density(double arg) const;




};


//
// inline members
//



inline
void
ExcitonProperties::set_statistics(TiberCad::Statistics statistics)
{
  _statistics = statistics;
}


inline
TiberCad::Statistics
ExcitonProperties::get_statistics(void) const
{
  return _statistics;
}

inline
double
ExcitonProperties::get_net_recombination_rate(void) const
{
  return net_recomb_rate;
}


inline
double
ExcitonProperties::get_nonradiative_recombination_rate(void)
{
  return 0;
}


inline
double
ExcitonProperties::get_radiative_recombination_rate(void)
{
  return 0;
}


inline
double
ExcitonProperties::get_dissociation_rate(void)
{
  return 0;
}


inline
double
ExcitonProperties::calculate_density(double arg) const
{

  const double arg_max = 150;
  const double arg_min = -100;

  double dens;
  if (arg < arg_max)
    dens = std::exp(arg);
  else
    dens = std::exp(arg_max);

  return dens;
}


inline
void
ExcitonProperties::calculate_density_and_derivative(double arg, double& density,
    double& derivative) const
{
  density = calculate_density(arg);
  derivative = density;
}


inline
const Elem*
ExcitonProperties::get_element(void) const
{
  return _elem;
}


inline
const Point&
ExcitonProperties::get_coordinates(void) const
{
  return *_coord;
}


inline
void
ExcitonProperties::set_coordinates(const Point& p)
{
  _coord = &p;
}


inline
void
ExcitonProperties::set_carrier_temperature(double T_x)
{
  exciton_vt = T_x;
}

inline
double
ExcitonProperties::get_lattice_temperature(void)
{
  return exciton_vt;
}

inline
double
ExcitonProperties::get_exciton_energy(void)
{
  return _energy;
}

inline
double
ExcitonProperties::get_carrier_temperature(void)
{
  return exciton_vt;
}



inline
void
ExcitonProperties::set_effective_potential(double eff_potential)
{
  eff_pot = eff_potential;
}


inline
void
ExcitonProperties::set_energy(double energy)
{
  _energy = energy;
}


inline
void
ExcitonProperties::set_density_of_states(double DOS)
{
  _DOS = DOS;
}


inline
void
ExcitonProperties::set_density(double x)
{
  density = x;
}


inline
void
ExcitonProperties::calculate_density(void)
{
  double kT = exciton_vt;

  double arg_x = - (_energy - eff_pot) / kT;

  calculate_density_and_derivative(arg_x, density, density_derivative);

  double Nx = _DOS;
  density *= Nx;
  density_derivative *= Nx / kT;
}


inline
void
ExcitonProperties::calculate_mobility(void)
{
  do_mobility();
}


inline
void
ExcitonProperties::calculate_net_recombination_rate(void)
{
  do_recombination();
}





#endif /* _EXCITONPROPERTIES_H_*/
