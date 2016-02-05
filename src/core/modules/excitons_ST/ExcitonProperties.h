// $Id: ExcitonProperties.h 4192 2015-12-10 11:11:18Z maufder $

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
     * \param mat the bulk material
     * \param options the options as given in the input file
     * \return a pointer to the newly created object
     */
    static ExcitonProperties* create(const std::string& name,
        const Material* mat, const ModelOptions& options = ModelOptions());

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

    //! Set the exciton energy \f$E_x\f$
    /*!
     * \f$E_x=E_g-R\f$ where \f$E_g\f$ and \f$R\f$ are the band gap and
     * the Exciton binding energy, respectively.
     */
    void set_energy(double energy);

    //! Set the exciton DOS
    void set_density_of_states(double DOS);

    //! Set the carrier density
    void set_density(double s, double t);


    //! Get the element we are currently working on
    const Elem* get_element(void) const;

    //! Get the coordinates of the point we are currently working on
    const Point& get_coordinates(void) const;


    //! Calculate the exciton net recombination rate
    void calculate_net_recombination_rate(void);

    //! Calculate the exciton diffusion
    void calculate_diffusion(void);


    //! Get the lattice temperature
    double get_lattice_temperature(void);

    //! Get the carrier temperature
    double get_carrier_temperature(void);

    //! Get the singlet exciton density
    double get_s_density(void) const
      { return s_density; };

    //! Get the triplet exciton density
    double get_t_density(void) const
      { return t_density; };


    //! Get the singlet recombination rate
    double get_s_net_recombination_rate(void) const;

    //! Get the triplet recombination rate
    double get_t_net_recombination_rate(void) const;

    //! Get the singlet recombination rate derivative
    /*!
     * Get \f$\frac{\partial R}{\partial\phi_x}\f$
     */
    double get_s_net_recombination_rate_derivative(void) const
        { return s_recombination_rate_derivative; };

    //! Get the singlet nonradiative recombination rate
    virtual double get_s_nonradiative_recombination_rate(void);

    //! Get the singlet radiative recombination rate
    virtual double get_s_radiative_recombination_rate(void);

    //! Get the singlet dissociation rate
    virtual double get_s_dissociation_rate(void);

    //! Get the ISC rate
    virtual double get_isc_rate(void);

    //! Get the ISC rate derivative;
    virtual double get_isc_rate_derivative(void);

    //! Get the triplet recombination rate derivative
    /*!
     * Get \f$\frac{\partial R}{\partial\phi_x}\f$
     */
    double get_t_net_recombination_rate_derivative(void) const
        { return t_recombination_rate_derivative; };

    //! Get the triplet nonradiative recombination rate
    virtual double get_t_nonradiative_recombination_rate(void);

    //! Get the triplet radiative recombination rate
    virtual double get_t_radiative_recombination_rate(void);

    //! Get the triplet dissociation rate
    virtual double get_t_dissociation_rate(void);

    //! Get the singlet mobility
    /*!
     * \return the singlet mobility
     */
    double get_s_diffusion(void) const
      { return s_diffusion; };

    //! Get the triplet mobility
    /*!
     * \return the triplet mobility
     */
    double get_t_diffusion(void) const
      { return t_diffusion; };

    virtual double get_s_generation_rate() {
      return 0;
    };

    virtual double get_t_generation_rate() {
      return 0;
    };

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

    virtual void do_diffusion(void) {};


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

    //! The singlets density
    double s_density;

    //! The triplets density
    double t_density;


    //! The singlet diffusion coefficient
    double s_diffusion;

    //! The singlet radiative recombination rate
    double s_rad_recomb_rate;

    //! The singlet nonradiative recombination rate
    double s_nonrad_recomb_rate;

    //! The singlet dissociation rate
    double s_dissoc_rate;

    //! The singlet net recombination rate
    double s_net_recomb_rate;

    //! The derivative of the singlet net recombination rate
    double s_recombination_rate_derivative;

    //! The triplet diffusion coefficient
    double t_diffusion;

    //! The triplet radiative recombination rate
    double t_rad_recomb_rate;

    //! The triplet nonradiative recombination rate
    double t_nonrad_recomb_rate;

    //! The triplet dissociation rate
    double t_dissoc_rate;

    //! The triplet net recombination rate
    double t_net_recomb_rate;

    //! The derivative of the triplet net recombination rate
    double t_recombination_rate_derivative;


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
ExcitonProperties::get_s_net_recombination_rate(void) const
{
  return s_net_recomb_rate;
}

inline
double
ExcitonProperties::get_s_nonradiative_recombination_rate(void)
{
  return 0;
}

inline
double
ExcitonProperties::get_s_radiative_recombination_rate(void)
{
  return 0;
}

inline
double
ExcitonProperties::get_s_dissociation_rate(void)
{
  return 0;
}

inline
double
ExcitonProperties::get_isc_rate(void)
{
  return 0;
}

inline
double
ExcitonProperties::get_isc_rate_derivative(void)
{
  return 0;
}

inline
double
ExcitonProperties::get_t_net_recombination_rate(void) const
{
  return t_net_recomb_rate;
}

inline
double
ExcitonProperties::get_t_nonradiative_recombination_rate(void)
{
  return 0;
}

inline
double
ExcitonProperties::get_t_radiative_recombination_rate(void)
{
  return 0;
}

inline
double
ExcitonProperties::get_t_dissociation_rate(void)
{
  return 0;
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
ExcitonProperties::set_density(double s, double t)
{
  s_density = s;
  t_density = t;
}


//inline
//void
//ExcitonProperties::calculate_density(void)
//{
  /*
  double kT = exciton_vt;

  double arg_x = - (_energy - eff_pot) / kT;

  calculate_density_and_derivative(arg_x, density, density_derivative);

  double Nx = _DOS;
  density *= Nx;
  density_derivative *= Nx / kT;
  */
//}


inline
void
ExcitonProperties::calculate_diffusion(void)
{
  do_diffusion();
}


inline
void
ExcitonProperties::calculate_net_recombination_rate(void)
{
  do_recombination();
}





#endif /* _EXCITONPROPERTIES_H_*/
